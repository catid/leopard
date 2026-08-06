/*
    Copyright (c) 2017 Christopher A. Taylor.  All rights reserved.
    See Leopard2Backend.h for the full BSD license notice.
*/

#include "Leopard2Backend.h"
#include "LeopardCommon.h"

#include <algorithm>
#include <cstring>
#include <immintrin.h>
#include <memory>
#include <new>

namespace leopard { namespace backend {

#if defined(LEO2_AVX512_VARIANT)
# if defined(LEO2_GFNI_VARIANT)
#  define LEO2_AVX_BACKEND_NAME "avx512vl-gfni"
# else
#  define LEO2_AVX_BACKEND_NAME "avx512vl"
# endif
# define LEO2_AVX_BACKEND_KIND LEO2_BACKEND_AVX512
# define LEO2_AVX_INITIALIZER InitializeAVX512
# define LEO2_AVX_TABLE_STATE TestGetAVX512TableState
#elif defined(LEO2_GFNI_MEMBER)
// Production GFNI member: its own backend identity, qualified by CPUID
// GFNI+AVX2, selected only by explicit request.  LEO2_GFNI_VARIANT without
// LEO2_GFNI_MEMBER remains the in-place evaluation configuration that reuses
// the AVX2 identity for A/B experiments.
# define LEO2_AVX_BACKEND_KIND LEO2_BACKEND_GFNI
# define LEO2_AVX_BACKEND_NAME "avx2-gfni"
# define LEO2_AVX_INITIALIZER InitializeGFNI
# define LEO2_AVX_TABLE_STATE TestGetGFNITableState
#else
# define LEO2_AVX_BACKEND_KIND LEO2_BACKEND_AVX2
# define LEO2_AVX_BACKEND_NAME "avx2"
# define LEO2_AVX_INITIALIZER InitializeAVX2
# define LEO2_AVX_TABLE_STATE TestGetAVX2TableState
#endif

// GFNI variant contract.
//
// A fixed multiplication by one field element is a GF(2)-linear map on the
// symbol bits, so each byte-wide component is an 8-by-8 bit matrix and
// VGF2P8AFFINEQB evaluates it in one instruction.  The measured operand order
// on this family is operand_bit = 8 * (7 - output_bit) + input_bit, which is
// the order the exhaustive experiment under experiments/leopard2/gfni_affine
// validated against all 65,536 GF8 products.
//
// GF8 deliberately reuses the nibble-table storage shape: each 16-byte row
// instead holds one affine matrix duplicated, and a 128-bit broadcast fills
// every 64-bit lane with that matrix.  GF16 uses four packed 64-bit matrices
// per multiplier and broadcasts them at use sites, reducing its table from
// 8 MiB to 2 MiB.  Scalar tails evaluate these same affine operands directly.
#if defined(LEO2_GFNI_VARIANT)
static LEO_FORCE_INLINE uint64_t GFNIAffineMatrixBit(
    unsigned output_bit, unsigned input_bit)
{
    return static_cast<uint64_t>(1) << (8 * (7 - output_bit) + input_bit);
}

#ifdef LEO_HAS_FF8
static void GFNIStoreMatrix(uint8_t row[16], uint64_t matrix)
{
    for (unsigned byte = 0; byte < 8; ++byte)
    {
        const uint8_t value = static_cast<uint8_t>(matrix >> (8 * byte));
        row[byte] = value;
        row[byte + 8] = value;
    }
}
#endif

static LEO_FORCE_INLINE uint8_t GFNIParity(uint8_t value)
{
    value = static_cast<uint8_t>(value ^ (value >> 4));
    value = static_cast<uint8_t>(value ^ (value >> 2));
    return static_cast<uint8_t>((value ^ (value >> 1)) & 1U);
}

// Scalar evaluation of one stored affine matrix.  Sub-vector tails use this
// instead of a second set of nibble tables, so the variant needs no extra
// storage beyond the affine operands themselves.  Operand byte 7 - output_bit
// holds that output bit's input mask, matching GFNIAffineMatrixBit.
static LEO_FORCE_INLINE uint8_t GFNIApplyMatrix(
    const uint8_t row[16], uint8_t value)
{
    uint8_t result = 0;
    for (unsigned output_bit = 0; output_bit < 8; ++output_bit)
        result = static_cast<uint8_t>(result |
            (GFNIParity(static_cast<uint8_t>(row[7 - output_bit] & value))
                << output_bit));
    return result;
}

// Packed-form twins for the GF16 affine tables, which store one 8-byte matrix
// per block instead of a duplicated 16-byte row.  Byte b of the integer is
// row byte b of the stored form above, so `set1_epi64x` fills the four 64-bit
// lanes exactly as the 16-byte broadcast did and VGF2P8AFFINEQB sees
// identical operands.
static LEO_FORCE_INLINE uint8_t GFNIApplyMatrix64(
    uint64_t matrix, uint8_t value)
{
    uint8_t result = 0;
    for (unsigned output_bit = 0; output_bit < 8; ++output_bit)
        result = static_cast<uint8_t>(result |
            (GFNIParity(static_cast<uint8_t>(
                static_cast<uint8_t>(matrix >> (8 * (7 - output_bit))) &
                value)) << output_bit));
    return result;
}

static LEO_FORCE_INLINE __m256i BroadcastAffine(uint64_t matrix)
{
    return _mm256_set1_epi64x(static_cast<long long>(matrix));
}
#endif

#ifdef LEO_HAS_FF8
struct FF8NibbleTable
{
    uint8_t low[16];
    uint8_t high[16];
};
static const FF8NibbleTable* FF8Tables = NULL;
#endif

#ifdef LEO_HAS_FF16
#if defined(LEO2_GFNI_VARIANT)
// Packed affine storage: the four 8x8 GF(2) matrix blocks of one fixed GF16
// multiplication, 32 bytes per logarithm (2 MiB total) instead of the nibble
// shape's 128 bytes (8 MiB).  Vector kernels broadcast each 64-bit matrix
// with vpbroadcastq, which fills the four 64-bit lanes exactly as the former
// 16-byte duplicated-row broadcast did, so VGF2P8AFFINEQB sees identical
// operands.  GFNI doc production requirement 2.
struct FF16AffineTable
{
    uint64_t block[4];
};
static const FF16AffineTable* FF16Tables = NULL;
typedef FF16AffineTable FF16Table;
#else
struct FF16NibbleTable
{
    uint8_t low[4][16];
    uint8_t high[4][16];
};
static const FF16NibbleTable* FF16Tables = NULL;
typedef FF16NibbleTable FF16Table;
#endif
#endif

#ifdef LEO_HAS_FF8
static uint8_t FF8Product(uint16_t log, uint8_t value)
{
    const FF8NibbleTable& table = FF8Tables[log];
#if defined(LEO2_GFNI_VARIANT)
    return GFNIApplyMatrix(table.low, value);
#else
    return static_cast<uint8_t>(
        table.low[value & 15U] ^ table.high[value >> 4]);
#endif
}
#endif

#if defined(LEO_HAS_FF8) && !defined(LEO2_AVX512_VARIANT) && \
    !defined(LEO2_GFNI_VARIANT)
static LEO_FORCE_INLINE __m256i AVX2AddMod255(
    __m256i a, __m256i b)
{
    const __m256i mask = _mm256_set1_epi16(255);
    __m256i sum = _mm256_add_epi16(a, b);
    sum = _mm256_add_epi16(sum, _mm256_srli_epi16(sum, 8));
    return _mm256_and_si256(sum, mask);
}

static LEO_FORCE_INLINE __m256i AVX2SubMod255(
    __m256i a, __m256i b)
{
    const __m256i modulus = _mm256_set1_epi16(255);
    const __m256i borrow = _mm256_cmpgt_epi16(b, a);
    return _mm256_add_epi16(
        _mm256_sub_epi16(a, b), _mm256_and_si256(borrow, modulus));
}

static LEO_FORCE_INLINE void AVX2WalshPair(
    __m256i a, __m256i b, __m256i& sum, __m256i& difference)
{
    sum = AVX2AddMod255(a, b);
    difference = AVX2SubMod255(a, b);
}

template<int Distance>
static LEO_FORCE_INLINE __m256i AVX2Walsh16Stage(__m256i value)
{
    __m256i swapped;
    if (Distance == 1)
    {
        swapped = _mm256_shufflelo_epi16(value, 0xb1);
        swapped = _mm256_shufflehi_epi16(swapped, 0xb1);
    }
    else if (Distance == 2)
    {
        swapped = _mm256_shufflelo_epi16(value, 0x4e);
        swapped = _mm256_shufflehi_epi16(swapped, 0x4e);
    }
    else if (Distance == 4)
        swapped = _mm256_shuffle_epi32(value, 0x4e);
    else
        swapped = _mm256_permute2x128_si256(value, value, 0x01);

    __m256i sum, difference;
    AVX2WalshPair(value, swapped, sum, difference);
    if (Distance == 1)
    {
        difference = _mm256_shufflelo_epi16(difference, 0xb1);
        difference = _mm256_shufflehi_epi16(difference, 0xb1);
        return _mm256_blend_epi16(sum, difference, 0xaa);
    }
    if (Distance == 2)
    {
        difference = _mm256_shufflelo_epi16(difference, 0x4e);
        difference = _mm256_shufflehi_epi16(difference, 0x4e);
        return _mm256_blend_epi16(sum, difference, 0xcc);
    }
    if (Distance == 4)
    {
        difference = _mm256_shuffle_epi32(difference, 0x4e);
        return _mm256_blend_epi16(sum, difference, 0xf0);
    }
    difference = _mm256_permute2x128_si256(
        difference, difference, 0x01);
    return _mm256_blend_epi32(sum, difference, 0xf0);
}

static LEO_FORCE_INLINE __m256i AVX2Walsh16(__m256i value)
{
    value = AVX2Walsh16Stage<1>(value);
    value = AVX2Walsh16Stage<2>(value);
    value = AVX2Walsh16Stage<4>(value);
    return AVX2Walsh16Stage<8>(value);
}

static LEO_FORCE_INLINE __m256i AVX2PackMod255(
    __m256i low, __m256i high)
{
    // vpackuswb operates independently in each 128-bit lane.  Restore the
    // original contiguous element order after packing.
    return _mm256_permute4x64_epi64(
        _mm256_packus_epi16(low, high), 0xd8);
}

static LEO_FORCE_INLINE void AVX2ExpandBytes(
    __m256i bytes, __m256i& low, __m256i& high)
{
    low = _mm256_cvtepu8_epi16(_mm256_castsi256_si128(bytes));
    high = _mm256_cvtepu8_epi16(_mm256_extracti128_si256(bytes, 1));
}

static void AVX2WalshTransform(uint8_t* data, uint32_t n)
{
    for (uint32_t offset = 0; offset < n; offset += 32)
    {
        const __m256i packed = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(data + offset));
        __m256i low, high;
        AVX2ExpandBytes(packed, low, high);
        low = AVX2Walsh16(low);
        high = AVX2Walsh16(high);
        AVX2WalshPair(low, high, low, high);
        _mm256_storeu_si256(reinterpret_cast<__m256i*>(data + offset),
            AVX2PackMod255(low, high));
    }

    for (uint32_t distance = 32; distance < n; distance <<= 1)
    {
        const uint32_t group_size = distance << 1;
        for (uint32_t group = 0; group < n; group += group_size)
        {
            for (uint32_t offset = 0; offset < distance; offset += 32)
            {
                uint8_t* const a_pointer = data + group + offset;
                uint8_t* const b_pointer = a_pointer + distance;
                __m256i a_low, a_high, b_low, b_high;
                AVX2ExpandBytes(_mm256_loadu_si256(
                    reinterpret_cast<const __m256i*>(a_pointer)),
                    a_low, a_high);
                AVX2ExpandBytes(_mm256_loadu_si256(
                    reinterpret_cast<const __m256i*>(b_pointer)),
                    b_low, b_high);
                __m256i sum_low, sum_high, difference_low, difference_high;
                AVX2WalshPair(
                    a_low, b_low, sum_low, difference_low);
                AVX2WalshPair(
                    a_high, b_high, sum_high, difference_high);
                _mm256_storeu_si256(
                    reinterpret_cast<__m256i*>(a_pointer),
                    AVX2PackMod255(sum_low, sum_high));
                _mm256_storeu_si256(
                    reinterpret_cast<__m256i*>(b_pointer),
                    AVX2PackMod255(difference_low, difference_high));
            }
        }
    }
}

static LEO_FORCE_INLINE __m256i AVX2MultiplyMod255(
    __m256i a, __m256i b)
{
    const __m256i mask = _mm256_set1_epi16(255);
    const __m256i cutoff = _mm256_set1_epi16(254);
    const __m256i product = _mm256_mullo_epi16(a, b);
    __m256i folded = _mm256_add_epi16(
        _mm256_and_si256(product, mask), _mm256_srli_epi16(product, 8));
    const __m256i reduce = _mm256_cmpgt_epi16(folded, cutoff);
    folded = _mm256_sub_epi16(folded, _mm256_and_si256(reduce, mask));
    return folded;
}

static void AVX2FF8WalshLocator(
    const uint8_t* erasures,
    const uint8_t* transformed_kernel,
    uint8_t* locator_logs,
    uint32_t n)
{
    const __m256i zero = _mm256_setzero_si256();
    const __m256i one = _mm256_set1_epi8(1);
    for (uint32_t offset = 0; offset < n; offset += 32)
    {
        __m256i erased = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(erasures + offset));
        erased = _mm256_andnot_si256(
            _mm256_cmpeq_epi8(erased, zero), one);
        _mm256_storeu_si256(
            reinterpret_cast<__m256i*>(locator_logs + offset), erased);
    }

    AVX2WalshTransform(locator_logs, n);
    for (uint32_t offset = 0; offset < n; offset += 32)
    {
        __m256i value_low, value_high, kernel_low, kernel_high;
        AVX2ExpandBytes(_mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(locator_logs + offset)),
            value_low, value_high);
        AVX2ExpandBytes(_mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(transformed_kernel + offset)),
            kernel_low, kernel_high);
        value_low = AVX2MultiplyMod255(value_low, kernel_low);
        value_high = AVX2MultiplyMod255(value_high, kernel_high);
        _mm256_storeu_si256(
            reinterpret_cast<__m256i*>(locator_logs + offset),
            AVX2PackMod255(value_low, value_high));
    }
    AVX2WalshTransform(locator_logs, n);
}
#endif

#if defined(LEO_HAS_FF8) || !defined(LEO2_GFNI_VARIANT)
static __m256i BroadcastTable(const uint8_t table[16])
{
    return _mm256_broadcastsi128_si256(_mm_loadu_si128(
        reinterpret_cast<const __m128i*>(table)));
}
#endif

#ifdef LEO_HAS_FF8
#if defined(LEO2_AVX512_VARIANT)
# define LEO2_AVX2_OPERATION_INLINE
#else
# define LEO2_AVX2_OPERATION_INLINE LEO_FORCE_INLINE
#endif

template<bool Add>
static LEO2_AVX2_OPERATION_INLINE void AVX2FF8Operation(
    void* destination,
    const void* source,
    uint16_t multiplier_log,
    uint64_t byte_count)
{
    uint8_t* output = static_cast<uint8_t*>(destination);
    const uint8_t* input = static_cast<const uint8_t*>(source);
    const FF8NibbleTable& table = FF8Tables[multiplier_log];
    const __m256i low_table = BroadcastTable(table.low);
#if !defined(LEO2_GFNI_VARIANT)
    const __m256i high_table = BroadcastTable(table.high);
    const __m256i nibble_mask = _mm256_set1_epi8(15);
#endif
#if defined(LEO2_GFNI_VARIANT)
    // The affine form leaves only one arithmetic operation per vector, so the
    // three loop-bookkeeping instructions would otherwise be a large share of
    // the issue slots.  Two independent vectors per iteration also give the
    // three-cycle affine latency somewhere to hide.  This matches Leopard1's
    // 64-byte mul_mem/muladd_mem stride.
    while (byte_count >= 64)
    {
        const __m256i data0 = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(input));
        const __m256i data1 = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(input + 32));
        __m256i product0 = _mm256_gf2p8affine_epi64_epi8(data0, low_table, 0);
        __m256i product1 = _mm256_gf2p8affine_epi64_epi8(data1, low_table, 0);
        if (Add)
        {
            product0 = _mm256_xor_si256(product0, _mm256_loadu_si256(
                reinterpret_cast<const __m256i*>(output)));
            product1 = _mm256_xor_si256(product1, _mm256_loadu_si256(
                reinterpret_cast<const __m256i*>(output + 32)));
        }
        _mm256_storeu_si256(reinterpret_cast<__m256i*>(output), product0);
        _mm256_storeu_si256(reinterpret_cast<__m256i*>(output + 32), product1);
        input += 64;
        output += 64;
        byte_count -= 64;
    }
#endif
    while (byte_count >= 32)
    {
        const __m256i data = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(input));
#if defined(LEO2_GFNI_VARIANT)
        __m256i product = _mm256_gf2p8affine_epi64_epi8(data, low_table, 0);
#else
        const __m256i low = _mm256_shuffle_epi8(
            low_table, _mm256_and_si256(data, nibble_mask));
        const __m256i high = _mm256_shuffle_epi8(high_table,
            _mm256_and_si256(_mm256_srli_epi64(data, 4), nibble_mask));
        __m256i product = _mm256_xor_si256(low, high);
#endif
        if (Add)
            product = _mm256_xor_si256(product, _mm256_loadu_si256(
                reinterpret_cast<const __m256i*>(output)));
        _mm256_storeu_si256(reinterpret_cast<__m256i*>(output), product);
        input += 32;
        output += 32;
        byte_count -= 32;
    }
    while (byte_count-- != 0)
    {
        const uint8_t product = FF8Product(multiplier_log, *input++);
        if (Add)
            *output++ ^= product;
        else
            *output++ = product;
    }
}

#undef LEO2_AVX2_OPERATION_INLINE


static void AVX2FF8Multiply(
    void* destination, const void* source,
    uint16_t multiplier_log, uint64_t byte_count)
{
    AVX2FF8Operation<false>(destination, source, multiplier_log, byte_count);
}

#if defined(LEO2_AVX512_VARIANT)
# define LEO2_AVX2_MULADD_NOINLINE
#elif defined(_MSC_VER)
# define LEO2_AVX2_MULADD_NOINLINE __declspec(noinline)
#elif defined(__GNUC__) || defined(__clang__)
# define LEO2_AVX2_MULADD_NOINLINE __attribute__((noinline))
#else
# define LEO2_AVX2_MULADD_NOINLINE
#endif

static LEO2_AVX2_MULADD_NOINLINE void AVX2FF8MultiplyAdd(
    void* destination, const void* source,
    uint16_t multiplier_log, uint64_t byte_count)
{
    AVX2FF8Operation<true>(destination, source, multiplier_log, byte_count);
}

#undef LEO2_AVX2_MULADD_NOINLINE

static __m256i AVX2FF8ProductVector(
    __m256i data,
    __m256i low_table,
    __m256i high_table)
{
#if defined(LEO2_GFNI_VARIANT)
    // low_table carries the broadcast affine matrix; high_table is unused.
    (void)high_table;
    return _mm256_gf2p8affine_epi64_epi8(data, low_table, 0);
#else
    const __m256i nibble_mask = _mm256_set1_epi8(15);
    const __m256i low = _mm256_shuffle_epi8(
        low_table, _mm256_and_si256(data, nibble_mask));
    const __m256i high = _mm256_shuffle_epi8(high_table,
        _mm256_and_si256(_mm256_srli_epi64(data, 4), nibble_mask));
    return _mm256_xor_si256(low, high);
#endif
}

static __m256i AVX2FF8ApplyWeight(__m256i data, uint16_t weight_log)
{
    static const uint16_t kModulus = 255;
    if (weight_log == 0 || weight_log == kModulus)
        return data;
    const FF8NibbleTable& table = FF8Tables[weight_log];
    return AVX2FF8ProductVector(
        data, BroadcastTable(table.low), BroadcastTable(table.high));
}

template<bool Inverse>
static void AVX2FF8Butterfly2(
    void* x_pointer,
    void* y_pointer,
    uint16_t multiplier_log,
    uint64_t byte_count)
{
    uint8_t* x = static_cast<uint8_t*>(x_pointer);
    uint8_t* y = static_cast<uint8_t*>(y_pointer);
    const FF8NibbleTable& table = FF8Tables[multiplier_log];
    const __m256i low_table = BroadcastTable(table.low);
    const __m256i high_table = BroadcastTable(table.high);
#if !defined(LEO2_AVX512_VARIANT)
    // Two independent vectors hide shuffle latency for small shards on the
    // explicitly qualified AVX2 backend.  Keep the measured 4 KiB cutoff: an
    // unconditional unroll regressed T8 at 64 KiB, while the mature
    // single-vector loop remains neutral above it.  The separately qualified
    // AVX-512VL variant retains its measured schedule.
    if (byte_count <= 4096)
    {
        while (byte_count >= 64)
        {
            __m256i x_value0 = _mm256_loadu_si256(
                reinterpret_cast<const __m256i*>(x));
            __m256i y_value0 = _mm256_loadu_si256(
                reinterpret_cast<const __m256i*>(y));
            __m256i x_value1 = _mm256_loadu_si256(
                reinterpret_cast<const __m256i*>(x + 32));
            __m256i y_value1 = _mm256_loadu_si256(
                reinterpret_cast<const __m256i*>(y + 32));
            if (Inverse)
            {
                y_value0 = _mm256_xor_si256(y_value0, x_value0);
                y_value1 = _mm256_xor_si256(y_value1, x_value1);
                x_value0 = _mm256_xor_si256(x_value0,
                    AVX2FF8ProductVector(y_value0, low_table, high_table));
                x_value1 = _mm256_xor_si256(x_value1,
                    AVX2FF8ProductVector(y_value1, low_table, high_table));
            }
            else
            {
                x_value0 = _mm256_xor_si256(x_value0,
                    AVX2FF8ProductVector(y_value0, low_table, high_table));
                x_value1 = _mm256_xor_si256(x_value1,
                    AVX2FF8ProductVector(y_value1, low_table, high_table));
                y_value0 = _mm256_xor_si256(y_value0, x_value0);
                y_value1 = _mm256_xor_si256(y_value1, x_value1);
            }
            _mm256_storeu_si256(reinterpret_cast<__m256i*>(x), x_value0);
            _mm256_storeu_si256(reinterpret_cast<__m256i*>(y), y_value0);
            _mm256_storeu_si256(
                reinterpret_cast<__m256i*>(x + 32), x_value1);
            _mm256_storeu_si256(
                reinterpret_cast<__m256i*>(y + 32), y_value1);
            x += 64;
            y += 64;
            byte_count -= 64;
        }
    }
#endif
    while (byte_count >= 32)
    {
        __m256i x_value = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(x));
        __m256i y_value = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(y));
        if (Inverse)
        {
            y_value = _mm256_xor_si256(y_value, x_value);
            x_value = _mm256_xor_si256(x_value,
                AVX2FF8ProductVector(y_value, low_table, high_table));
        }
        else
        {
            x_value = _mm256_xor_si256(x_value,
                AVX2FF8ProductVector(y_value, low_table, high_table));
            y_value = _mm256_xor_si256(y_value, x_value);
        }
        _mm256_storeu_si256(reinterpret_cast<__m256i*>(x), x_value);
        _mm256_storeu_si256(reinterpret_cast<__m256i*>(y), y_value);
        x += 32;
        y += 32;
        byte_count -= 32;
    }
    while (byte_count-- != 0)
    {
        if (Inverse)
        {
            *y ^= *x;
            *x ^= FF8Product(multiplier_log, *y);
        }
        else
        {
            *x ^= FF8Product(multiplier_log, *y);
            *y ^= *x;
        }
        ++x;
        ++y;
    }
}

static void AVX2FF8IFFTButterfly2(
    void* x, void* y, uint16_t multiplier_log, uint64_t byte_count)
{
    AVX2FF8Butterfly2<true>(x, y, multiplier_log, byte_count);
}

static void AVX2FF8FFTButterfly2(
    void* x, void* y, uint16_t multiplier_log, uint64_t byte_count)
{
    AVX2FF8Butterfly2<false>(x, y, multiplier_log, byte_count);
}

static void AVX2FF8FFTButterfly2Out(
    const void* x_input_pointer,
    const void* y_input_pointer,
    void* x_output_pointer,
    void* y_output_pointer,
    uint16_t multiplier_log,
    uint64_t byte_count)
{
    static const uint16_t kZeroSkew = 255;
    const uint8_t* x_input = static_cast<const uint8_t*>(x_input_pointer);
    const uint8_t* y_input = static_cast<const uint8_t*>(y_input_pointer);
    uint8_t* x_output = static_cast<uint8_t*>(x_output_pointer);
    uint8_t* y_output = static_cast<uint8_t*>(y_output_pointer);
    __m256i low_table = _mm256_setzero_si256();
    __m256i high_table = _mm256_setzero_si256();
    if (multiplier_log != kZeroSkew)
    {
        low_table = BroadcastTable(FF8Tables[multiplier_log].low);
        high_table = BroadcastTable(FF8Tables[multiplier_log].high);
    }
    while (byte_count >= 32)
    {
        __m256i x_value = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(x_input));
        __m256i y_value = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(y_input));
        if (multiplier_log != kZeroSkew)
            x_value = _mm256_xor_si256(x_value,
                AVX2FF8ProductVector(y_value, low_table, high_table));
        y_value = _mm256_xor_si256(y_value, x_value);
        _mm256_storeu_si256(reinterpret_cast<__m256i*>(x_output), x_value);
        _mm256_storeu_si256(reinterpret_cast<__m256i*>(y_output), y_value);
        x_input += 32;
        y_input += 32;
        x_output += 32;
        y_output += 32;
        byte_count -= 32;
    }
    while (byte_count-- != 0)
    {
        uint8_t x_value = *x_input++;
        uint8_t y_value = *y_input++;
        if (multiplier_log != kZeroSkew)
            x_value ^= FF8Product(multiplier_log, y_value);
        y_value ^= x_value;
        *x_output++ = x_value;
        *y_output++ = y_value;
    }
}

#if !defined(LEO2_AVX512_VARIANT)
static void AVX2FF8IFFTButterfly2Out(
    const void* x_input_pointer,
    const void* y_input_pointer,
    void* x_output_pointer,
    void* y_output_pointer,
    uint16_t multiplier_log,
    uint64_t byte_count)
{
    static const uint16_t kZeroSkew = 255;
    const uint8_t* x_input = static_cast<const uint8_t*>(x_input_pointer);
    const uint8_t* y_input = static_cast<const uint8_t*>(y_input_pointer);
    uint8_t* x_output = static_cast<uint8_t*>(x_output_pointer);
    uint8_t* y_output = static_cast<uint8_t*>(y_output_pointer);
    __m256i low_table = _mm256_setzero_si256();
    __m256i high_table = _mm256_setzero_si256();
    if (multiplier_log != kZeroSkew)
    {
        low_table = BroadcastTable(FF8Tables[multiplier_log].low);
        high_table = BroadcastTable(FF8Tables[multiplier_log].high);
    }
    while (byte_count >= 32)
    {
        const __m256i x_original = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(x_input));
        const __m256i y_value = _mm256_xor_si256(x_original,
            _mm256_loadu_si256(reinterpret_cast<const __m256i*>(y_input)));
        __m256i x_value = x_original;
        if (multiplier_log != kZeroSkew)
            x_value = _mm256_xor_si256(x_value,
                AVX2FF8ProductVector(y_value, low_table, high_table));
        _mm256_storeu_si256(reinterpret_cast<__m256i*>(x_output), x_value);
        _mm256_storeu_si256(reinterpret_cast<__m256i*>(y_output), y_value);
        x_input += 32;
        y_input += 32;
        x_output += 32;
        y_output += 32;
        byte_count -= 32;
    }
    while (byte_count-- != 0)
    {
        const uint8_t y_value = static_cast<uint8_t>(*y_input ^ *x_input);
        uint8_t x_value = *x_input;
        if (multiplier_log != kZeroSkew)
            x_value ^= FF8Product(multiplier_log, y_value);
        *x_output++ = x_value;
        *y_output++ = y_value;
        ++x_input;
        ++y_input;
    }
}
#endif

static void AVX2FF8IFFTButterfly2Xor(
    const void* x_input_pointer,
    const void* y_input_pointer,
    void* x_output_pointer,
    void* y_output_pointer,
    uint16_t multiplier_log,
    uint64_t byte_count)
{
    const uint8_t* x_input = static_cast<const uint8_t*>(x_input_pointer);
    const uint8_t* y_input = static_cast<const uint8_t*>(y_input_pointer);
    uint8_t* x_output = static_cast<uint8_t*>(x_output_pointer);
    uint8_t* y_output = static_cast<uint8_t*>(y_output_pointer);
    const FF8NibbleTable& table = FF8Tables[multiplier_log];
    const __m256i low_table = BroadcastTable(table.low);
    const __m256i high_table = BroadcastTable(table.high);
    while (byte_count >= 32)
    {
        const __m256i x_original = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(x_input));
        const __m256i y_value = _mm256_xor_si256(x_original,
            _mm256_loadu_si256(reinterpret_cast<const __m256i*>(y_input)));
        const __m256i x_value = _mm256_xor_si256(x_original,
            AVX2FF8ProductVector(y_value, low_table, high_table));
        const __m256i x_result = _mm256_xor_si256(x_value,
            _mm256_loadu_si256(reinterpret_cast<const __m256i*>(x_output)));
        const __m256i y_result = _mm256_xor_si256(y_value,
            _mm256_loadu_si256(reinterpret_cast<const __m256i*>(y_output)));
        _mm256_storeu_si256(reinterpret_cast<__m256i*>(x_output), x_result);
        _mm256_storeu_si256(reinterpret_cast<__m256i*>(y_output), y_result);
        x_input += 32;
        y_input += 32;
        x_output += 32;
        y_output += 32;
        byte_count -= 32;
    }
    while (byte_count-- != 0)
    {
        const uint8_t y_value = static_cast<uint8_t>(*y_input ^ *x_input);
        const uint8_t x_value = static_cast<uint8_t>(
            *x_input ^ FF8Product(multiplier_log, y_value));
        *x_output ^= x_value;
        *y_output ^= y_value;
        ++x_input;
        ++y_input;
        ++x_output;
        ++y_output;
    }
}

#endif // LEO_HAS_FF8

#ifdef LEO_HAS_FF16
static uint16_t FF16Product(uint16_t log, uint16_t value)
{
    const FF16Table& table = FF16Tables[log];
#if defined(LEO2_GFNI_VARIANT)
    const uint8_t input_low = static_cast<uint8_t>(value);
    const uint8_t input_high = static_cast<uint8_t>(value >> 8);
    const uint8_t product_low = static_cast<uint8_t>(
        GFNIApplyMatrix64(table.block[0], input_low) ^
        GFNIApplyMatrix64(table.block[1], input_high));
    const uint8_t product_high = static_cast<uint8_t>(
        GFNIApplyMatrix64(table.block[2], input_low) ^
        GFNIApplyMatrix64(table.block[3], input_high));
    return static_cast<uint16_t>(
        product_low | (static_cast<unsigned>(product_high) << 8));
#else
    const unsigned n0 = value & 15U;
    const unsigned n1 = (value >> 4) & 15U;
    const unsigned n2 = (value >> 8) & 15U;
    const unsigned n3 = value >> 12;
    const uint8_t low = static_cast<uint8_t>(
        table.low[0][n0] ^ table.low[1][n1] ^
        table.low[2][n2] ^ table.low[3][n3]);
    const uint8_t high = static_cast<uint8_t>(
        table.high[0][n0] ^ table.high[1][n1] ^
        table.high[2][n2] ^ table.high[3][n3]);
    return static_cast<uint16_t>(low | (static_cast<unsigned>(high) << 8));
#endif
}

static void AVX2FF16ProductVectors(
    __m256i low_data,
    __m256i high_data,
    const __m256i low_tables[4],
    const __m256i high_tables[4],
    __m256i& product_low,
    __m256i& product_high)
{
#if defined(LEO2_GFNI_VARIANT)
    // low_tables carries the four broadcast affine blocks of the 16-by-16 bit
    // multiplication matrix, ordered (low<-low, low<-high, high<-low,
    // high<-high).  high_tables is unused.
    (void)high_tables;
    product_low = _mm256_xor_si256(
        _mm256_gf2p8affine_epi64_epi8(low_data, low_tables[0], 0),
        _mm256_gf2p8affine_epi64_epi8(high_data, low_tables[1], 0));
    product_high = _mm256_xor_si256(
        _mm256_gf2p8affine_epi64_epi8(low_data, low_tables[2], 0),
        _mm256_gf2p8affine_epi64_epi8(high_data, low_tables[3], 0));
#else
    const __m256i mask = _mm256_set1_epi8(15);
    const __m256i nibbles[4] = {
        _mm256_and_si256(low_data, mask),
        _mm256_and_si256(_mm256_srli_epi64(low_data, 4), mask),
        _mm256_and_si256(high_data, mask),
        _mm256_and_si256(_mm256_srli_epi64(high_data, 4), mask)
    };
    product_low = _mm256_setzero_si256();
    product_high = _mm256_setzero_si256();
    for (unsigned i = 0; i < 4; ++i)
    {
        product_low = _mm256_xor_si256(product_low,
            _mm256_shuffle_epi8(low_tables[i], nibbles[i]));
        product_high = _mm256_xor_si256(product_high,
            _mm256_shuffle_epi8(high_tables[i], nibbles[i]));
    }
#endif
}

template<bool Add>
static void AVX2FF16Operation(
    void* destination,
    const void* source,
    uint16_t multiplier_log,
    uint64_t byte_count)
{
    uint8_t* output = static_cast<uint8_t*>(destination);
    const uint8_t* input = static_cast<const uint8_t*>(source);
    const FF16Table& table = FF16Tables[multiplier_log];
#if defined(LEO2_GFNI_VARIANT)
    const __m256i low_tables[4] = {
        BroadcastAffine(table.block[0]), BroadcastAffine(table.block[1]),
        BroadcastAffine(table.block[2]), BroadcastAffine(table.block[3])
    };
#else
    const __m256i low_tables[4] = {
        BroadcastTable(table.low[0]), BroadcastTable(table.low[1]),
        BroadcastTable(table.low[2]), BroadcastTable(table.low[3])
    };
#endif
#if !defined(LEO2_GFNI_VARIANT)
    const __m256i high_tables[4] = {
        BroadcastTable(table.high[0]), BroadcastTable(table.high[1]),
        BroadcastTable(table.high[2]), BroadcastTable(table.high[3])
    };
    const __m256i mask = _mm256_set1_epi8(15);
#endif
    uint64_t offset = 0;
    while (byte_count - offset >= 64)
    {
        const __m256i low_data = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(input + offset));
        const __m256i high_data = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(input + offset + 32));
#if defined(LEO2_GFNI_VARIANT)
        __m256i product_low = _mm256_xor_si256(
            _mm256_gf2p8affine_epi64_epi8(low_data, low_tables[0], 0),
            _mm256_gf2p8affine_epi64_epi8(high_data, low_tables[1], 0));
        __m256i product_high = _mm256_xor_si256(
            _mm256_gf2p8affine_epi64_epi8(low_data, low_tables[2], 0),
            _mm256_gf2p8affine_epi64_epi8(high_data, low_tables[3], 0));
#else
        const __m256i nibbles[4] = {
            _mm256_and_si256(low_data, mask),
            _mm256_and_si256(_mm256_srli_epi64(low_data, 4), mask),
            _mm256_and_si256(high_data, mask),
            _mm256_and_si256(_mm256_srli_epi64(high_data, 4), mask)
        };
        __m256i product_low = _mm256_setzero_si256();
        __m256i product_high = _mm256_setzero_si256();
        for (unsigned i = 0; i < 4; ++i)
        {
            product_low = _mm256_xor_si256(product_low,
                _mm256_shuffle_epi8(low_tables[i], nibbles[i]));
            product_high = _mm256_xor_si256(product_high,
                _mm256_shuffle_epi8(high_tables[i], nibbles[i]));
        }
#endif
        if (Add)
        {
            product_low = _mm256_xor_si256(product_low,
                _mm256_loadu_si256(reinterpret_cast<const __m256i*>(
                    output + offset)));
            product_high = _mm256_xor_si256(product_high,
                _mm256_loadu_si256(reinterpret_cast<const __m256i*>(
                    output + offset + 32)));
        }
        _mm256_storeu_si256(
            reinterpret_cast<__m256i*>(output + offset), product_low);
        _mm256_storeu_si256(
            reinterpret_cast<__m256i*>(output + offset + 32), product_high);
        offset += 64;
    }

    const uint64_t symbols = (byte_count - offset) / 2;
    for (uint64_t i = 0; i < symbols; ++i)
    {
        const uint16_t value = static_cast<uint16_t>(input[offset + i] |
            (static_cast<unsigned>(input[offset + symbols + i]) << 8));
        const uint16_t product = FF16Product(multiplier_log, value);
        if (Add)
        {
            output[offset + i] ^= static_cast<uint8_t>(product);
            output[offset + symbols + i] ^=
                static_cast<uint8_t>(product >> 8);
        }
        else
        {
            output[offset + i] = static_cast<uint8_t>(product);
            output[offset + symbols + i] = static_cast<uint8_t>(product >> 8);
        }
    }
}

static void AVX2FF16Multiply(
    void* destination, const void* source,
    uint16_t multiplier_log, uint64_t byte_count)
{
    AVX2FF16Operation<false>(destination, source, multiplier_log, byte_count);
}

static void AVX2FF16MultiplyAdd(
    void* destination, const void* source,
    uint16_t multiplier_log, uint64_t byte_count)
{
    AVX2FF16Operation<true>(destination, source, multiplier_log, byte_count);
}

template<bool Inverse>
static LEO_FORCE_INLINE void AVX2FF16Butterfly2Prepared(
    void* x_pointer,
    void* y_pointer,
    uint16_t multiplier_log,
    const __m256i low_tables[4],
    const __m256i high_tables[4],
    uint64_t byte_count)
{
    uint8_t* x = static_cast<uint8_t*>(x_pointer);
    uint8_t* y = static_cast<uint8_t*>(y_pointer);
    uint64_t offset = 0;
    while (byte_count - offset >= 64)
    {
        __m256i x_low = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(x + offset));
        __m256i x_high = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(x + offset + 32));
        __m256i y_low = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(y + offset));
        __m256i y_high = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(y + offset + 32));
        if (Inverse)
        {
            y_low = _mm256_xor_si256(y_low, x_low);
            y_high = _mm256_xor_si256(y_high, x_high);
        }
        __m256i product_low;
        __m256i product_high;
        AVX2FF16ProductVectors(y_low, y_high, low_tables, high_tables,
            product_low, product_high);
        x_low = _mm256_xor_si256(x_low, product_low);
        x_high = _mm256_xor_si256(x_high, product_high);
        if (!Inverse)
        {
            y_low = _mm256_xor_si256(y_low, x_low);
            y_high = _mm256_xor_si256(y_high, x_high);
        }
        _mm256_storeu_si256(reinterpret_cast<__m256i*>(x + offset), x_low);
        _mm256_storeu_si256(
            reinterpret_cast<__m256i*>(x + offset + 32), x_high);
        _mm256_storeu_si256(reinterpret_cast<__m256i*>(y + offset), y_low);
        _mm256_storeu_si256(
            reinterpret_cast<__m256i*>(y + offset + 32), y_high);
        offset += 64;
    }
    const uint64_t symbols = (byte_count - offset) / 2;
    for (uint64_t i = 0; i < symbols; ++i)
    {
        uint16_t x_value = static_cast<uint16_t>(x[offset + i] |
            (static_cast<unsigned>(x[offset + symbols + i]) << 8));
        uint16_t y_value = static_cast<uint16_t>(y[offset + i] |
            (static_cast<unsigned>(y[offset + symbols + i]) << 8));
        if (Inverse)
        {
            y_value ^= x_value;
            x_value ^= FF16Product(multiplier_log, y_value);
        }
        else
        {
            x_value ^= FF16Product(multiplier_log, y_value);
            y_value ^= x_value;
        }
        x[offset + i] = static_cast<uint8_t>(x_value);
        x[offset + symbols + i] = static_cast<uint8_t>(x_value >> 8);
        y[offset + i] = static_cast<uint8_t>(y_value);
        y[offset + symbols + i] = static_cast<uint8_t>(y_value >> 8);
    }
}

template<bool Inverse>
static void AVX2FF16Butterfly2(
    void* x_pointer,
    void* y_pointer,
    uint16_t multiplier_log,
    uint64_t byte_count)
{
    const FF16Table& table = FF16Tables[multiplier_log];
#if defined(LEO2_GFNI_VARIANT)
    const __m256i low_tables[4] = {
        BroadcastAffine(table.block[0]), BroadcastAffine(table.block[1]),
        BroadcastAffine(table.block[2]), BroadcastAffine(table.block[3])
    };
    // ProductVectors ignores its high_tables operand in this variant.
    const __m256i* const high_tables = low_tables;
#else
    const __m256i low_tables[4] = {
        BroadcastTable(table.low[0]), BroadcastTable(table.low[1]),
        BroadcastTable(table.low[2]), BroadcastTable(table.low[3])
    };
    const __m256i high_tables[4] = {
        BroadcastTable(table.high[0]), BroadcastTable(table.high[1]),
        BroadcastTable(table.high[2]), BroadcastTable(table.high[3])
    };
#endif
    AVX2FF16Butterfly2Prepared<Inverse>(
        x_pointer, y_pointer, multiplier_log,
        low_tables, high_tables, byte_count);
}

static void AVX2FF16IFFTButterfly2(
    void* x, void* y, uint16_t multiplier_log, uint64_t byte_count)
{
    AVX2FF16Butterfly2<true>(x, y, multiplier_log, byte_count);
}

static void AVX2FF16FFTButterfly2(
    void* x, void* y, uint16_t multiplier_log, uint64_t byte_count)
{
    AVX2FF16Butterfly2<false>(x, y, multiplier_log, byte_count);
}

static void AVX2FF16FFTButterfly2Out(
    const void* x_input_pointer,
    const void* y_input_pointer,
    void* x_output_pointer,
    void* y_output_pointer,
    uint16_t multiplier_log,
    uint64_t byte_count)
{
    static const uint16_t kZeroSkew = 65535;
    const uint8_t* x_input = static_cast<const uint8_t*>(x_input_pointer);
    const uint8_t* y_input = static_cast<const uint8_t*>(y_input_pointer);
    uint8_t* x_output = static_cast<uint8_t*>(x_output_pointer);
    uint8_t* y_output = static_cast<uint8_t*>(y_output_pointer);
    __m256i low_tables[4];
    __m256i high_tables[4];
    if (multiplier_log != kZeroSkew)
    {
        const FF16Table& table = FF16Tables[multiplier_log];
        for (unsigned i = 0; i < 4; ++i)
        {
#if defined(LEO2_GFNI_VARIANT)
            low_tables[i] = BroadcastAffine(table.block[i]);
            high_tables[i] = low_tables[i];
#else
            low_tables[i] = BroadcastTable(table.low[i]);
            high_tables[i] = BroadcastTable(table.high[i]);
#endif
        }
    }
    uint64_t offset = 0;
    while (byte_count - offset >= 64)
    {
        __m256i x_low = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(x_input + offset));
        __m256i x_high = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(x_input + offset + 32));
        __m256i y_low = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(y_input + offset));
        __m256i y_high = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(y_input + offset + 32));
        if (multiplier_log != kZeroSkew)
        {
            __m256i product_low;
            __m256i product_high;
            AVX2FF16ProductVectors(y_low, y_high, low_tables, high_tables,
                product_low, product_high);
            x_low = _mm256_xor_si256(x_low, product_low);
            x_high = _mm256_xor_si256(x_high, product_high);
        }
        y_low = _mm256_xor_si256(y_low, x_low);
        y_high = _mm256_xor_si256(y_high, x_high);
        _mm256_storeu_si256(
            reinterpret_cast<__m256i*>(x_output + offset), x_low);
        _mm256_storeu_si256(
            reinterpret_cast<__m256i*>(x_output + offset + 32), x_high);
        _mm256_storeu_si256(
            reinterpret_cast<__m256i*>(y_output + offset), y_low);
        _mm256_storeu_si256(
            reinterpret_cast<__m256i*>(y_output + offset + 32), y_high);
        offset += 64;
    }
    const uint64_t symbols = (byte_count - offset) / 2;
    for (uint64_t i = 0; i < symbols; ++i)
    {
        uint16_t x_value = static_cast<uint16_t>(x_input[offset + i] |
            (static_cast<unsigned>(x_input[offset + symbols + i]) << 8));
        uint16_t y_value = static_cast<uint16_t>(y_input[offset + i] |
            (static_cast<unsigned>(y_input[offset + symbols + i]) << 8));
        if (multiplier_log != kZeroSkew)
            x_value ^= FF16Product(multiplier_log, y_value);
        y_value ^= x_value;
        x_output[offset + i] = static_cast<uint8_t>(x_value);
        x_output[offset + symbols + i] = static_cast<uint8_t>(x_value >> 8);
        y_output[offset + i] = static_cast<uint8_t>(y_value);
        y_output[offset + symbols + i] = static_cast<uint8_t>(y_value >> 8);
    }
}

static void AVX2FF16IFFTButterfly2Xor(
    const void* x_input_pointer,
    const void* y_input_pointer,
    void* x_output_pointer,
    void* y_output_pointer,
    uint16_t multiplier_log,
    uint64_t byte_count)
{
    const uint8_t* x_input = static_cast<const uint8_t*>(x_input_pointer);
    const uint8_t* y_input = static_cast<const uint8_t*>(y_input_pointer);
    uint8_t* x_output = static_cast<uint8_t*>(x_output_pointer);
    uint8_t* y_output = static_cast<uint8_t*>(y_output_pointer);
    const FF16Table& table = FF16Tables[multiplier_log];
#if defined(LEO2_GFNI_VARIANT)
    const __m256i low_tables[4] = {
        BroadcastAffine(table.block[0]), BroadcastAffine(table.block[1]),
        BroadcastAffine(table.block[2]), BroadcastAffine(table.block[3])
    };
    // ProductVectors ignores its high_tables operand in this variant.
    const __m256i* const high_tables = low_tables;
#else
    const __m256i low_tables[4] = {
        BroadcastTable(table.low[0]), BroadcastTable(table.low[1]),
        BroadcastTable(table.low[2]), BroadcastTable(table.low[3])
    };
    const __m256i high_tables[4] = {
        BroadcastTable(table.high[0]), BroadcastTable(table.high[1]),
        BroadcastTable(table.high[2]), BroadcastTable(table.high[3])
    };
#endif
    uint64_t offset = 0;
    while (byte_count - offset >= 64)
    {
        __m256i x_low = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(x_input + offset));
        __m256i x_high = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(x_input + offset + 32));
        __m256i y_low = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(y_input + offset));
        __m256i y_high = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(y_input + offset + 32));
        y_low = _mm256_xor_si256(y_low, x_low);
        y_high = _mm256_xor_si256(y_high, x_high);
        __m256i product_low;
        __m256i product_high;
        AVX2FF16ProductVectors(y_low, y_high, low_tables, high_tables,
            product_low, product_high);
        x_low = _mm256_xor_si256(x_low, product_low);
        x_high = _mm256_xor_si256(x_high, product_high);
        x_low = _mm256_xor_si256(x_low, _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(x_output + offset)));
        x_high = _mm256_xor_si256(x_high, _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(x_output + offset + 32)));
        y_low = _mm256_xor_si256(y_low, _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(y_output + offset)));
        y_high = _mm256_xor_si256(y_high, _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(y_output + offset + 32)));
        _mm256_storeu_si256(
            reinterpret_cast<__m256i*>(x_output + offset), x_low);
        _mm256_storeu_si256(
            reinterpret_cast<__m256i*>(x_output + offset + 32), x_high);
        _mm256_storeu_si256(
            reinterpret_cast<__m256i*>(y_output + offset), y_low);
        _mm256_storeu_si256(
            reinterpret_cast<__m256i*>(y_output + offset + 32), y_high);
        offset += 64;
    }
    const uint64_t symbols = (byte_count - offset) / 2;
    for (uint64_t i = 0; i < symbols; ++i)
    {
        uint16_t x_value = static_cast<uint16_t>(x_input[offset + i] |
            (static_cast<unsigned>(x_input[offset + symbols + i]) << 8));
        uint16_t y_value = static_cast<uint16_t>(y_input[offset + i] |
            (static_cast<unsigned>(y_input[offset + symbols + i]) << 8));
        y_value ^= x_value;
        x_value ^= FF16Product(multiplier_log, y_value);
        x_output[offset + i] ^= static_cast<uint8_t>(x_value);
        x_output[offset + symbols + i] ^=
            static_cast<uint8_t>(x_value >> 8);
        y_output[offset + i] ^= static_cast<uint8_t>(y_value);
        y_output[offset + symbols + i] ^=
            static_cast<uint8_t>(y_value >> 8);
    }
}

#endif // LEO_HAS_FF16

static void AVX2XorMemory(
    void* destination, const void* source, uint64_t byte_count)
{
    uint8_t* output = static_cast<uint8_t*>(destination);
    const uint8_t* input = static_cast<const uint8_t*>(source);
    while (byte_count >= 32)
    {
        const __m256i result = _mm256_xor_si256(
            _mm256_loadu_si256(reinterpret_cast<const __m256i*>(output)),
            _mm256_loadu_si256(reinterpret_cast<const __m256i*>(input)));
        _mm256_storeu_si256(reinterpret_cast<__m256i*>(output), result);
        output += 32;
        input += 32;
        byte_count -= 32;
    }
    while (byte_count-- != 0)
        *output++ ^= *input++;
}

static void AVX2XorMemory2To1(
    void* destination,
    const void* source0,
    const void* source1,
    uint64_t byte_count)
{
    uint8_t* output = static_cast<uint8_t*>(destination);
    const uint8_t* input0 = static_cast<const uint8_t*>(source0);
    const uint8_t* input1 = static_cast<const uint8_t*>(source1);
    while (byte_count >= 128)
    {
        for (unsigned offset = 0; offset < 128; offset += 32)
        {
            __m256i result = _mm256_xor_si256(
                _mm256_loadu_si256(reinterpret_cast<const __m256i*>(
                    output + offset)),
                _mm256_loadu_si256(reinterpret_cast<const __m256i*>(
                    input0 + offset)));
            result = _mm256_xor_si256(result,
                _mm256_loadu_si256(reinterpret_cast<const __m256i*>(
                    input1 + offset)));
            _mm256_storeu_si256(
                reinterpret_cast<__m256i*>(output + offset), result);
        }
        output += 128;
        input0 += 128;
        input1 += 128;
        byte_count -= 128;
    }
    while (byte_count >= 32)
    {
        __m256i result = _mm256_xor_si256(
            _mm256_loadu_si256(reinterpret_cast<const __m256i*>(output)),
            _mm256_loadu_si256(reinterpret_cast<const __m256i*>(input0)));
        result = _mm256_xor_si256(result,
            _mm256_loadu_si256(reinterpret_cast<const __m256i*>(input1)));
        _mm256_storeu_si256(reinterpret_cast<__m256i*>(output), result);
        output += 32;
        input0 += 32;
        input1 += 32;
        byte_count -= 32;
    }
    while (byte_count-- != 0)
        *output++ ^= *input0++ ^ *input1++;
}

#if !defined(LEO2_AVX512_VARIANT)
template<unsigned SourceCount>
static void AVX2XorMemorySourceGroup(
    void* destination,
    const void* initial_source,
    const void* const* sources,
    uint64_t byte_count)
{
    uint8_t* output = static_cast<uint8_t*>(destination);
    const uint8_t* initial = static_cast<const uint8_t*>(initial_source);
    const uint8_t* inputs[SourceCount];
    for (unsigned lane = 0; lane < SourceCount; ++lane)
        inputs[lane] = static_cast<const uint8_t*>(sources[lane]);
    uint64_t offset = 0;
    while (byte_count - offset >= 32)
    {
        __m256i result = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(initial + offset));
        for (unsigned lane = 0; lane < SourceCount; ++lane)
            result = _mm256_xor_si256(result, _mm256_loadu_si256(
                reinterpret_cast<const __m256i*>(inputs[lane] + offset)));
        _mm256_storeu_si256(
            reinterpret_cast<__m256i*>(output + offset), result);
        offset += 32;
    }

    /* The ragged suffix is at most 31 bytes.  Volatile prevents the compiler
       from cloning large overlap-checking vector loops for this cold tail;
       aligned shard traffic retains the fully unrolled AVX2 loop above. */
    volatile uint8_t* const tail_output = output;
    const volatile uint8_t* const tail_initial = initial;
    while (offset < byte_count)
    {
        uint8_t result = tail_initial[offset];
        for (unsigned lane = 0; lane < SourceCount; ++lane)
            result ^= static_cast<const volatile uint8_t*>(
                inputs[lane])[offset];
        tail_output[offset++] = result;
    }
}

template<unsigned SourceCount>
static void AVX2XorMemoryDenseGroup(
    void* destination,
    const void* const* sources,
    uint64_t byte_count)
{
    uint8_t* output = static_cast<uint8_t*>(destination);
    const uint8_t* inputs[SourceCount];
    for (unsigned lane = 0; lane < SourceCount; ++lane)
        inputs[lane] = static_cast<const uint8_t*>(sources[lane]);
    uint64_t offset = 0;
    while (byte_count - offset >= 128)
    {
        for (unsigned block = 0; block < 128; block += 32)
        {
            __m256i result = _mm256_loadu_si256(
                reinterpret_cast<const __m256i*>(inputs[0] + offset + block));
            for (unsigned lane = 1; lane < SourceCount; ++lane)
                result = _mm256_xor_si256(result, _mm256_loadu_si256(
                    reinterpret_cast<const __m256i*>(
                        inputs[lane] + offset + block)));
            _mm256_storeu_si256(
                reinterpret_cast<__m256i*>(output + offset + block), result);
        }
        offset += 128;
    }
    while (byte_count - offset >= 32)
    {
        __m256i result = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(inputs[0] + offset));
        for (unsigned lane = 1; lane < SourceCount; ++lane)
            result = _mm256_xor_si256(result, _mm256_loadu_si256(
                reinterpret_cast<const __m256i*>(inputs[lane] + offset)));
        _mm256_storeu_si256(
            reinterpret_cast<__m256i*>(output + offset), result);
        offset += 32;
    }
    while (offset < byte_count)
    {
        uint8_t result = inputs[0][offset];
        for (unsigned lane = 1; lane < SourceCount; ++lane)
            result ^= inputs[lane][offset];
        output[offset++] = result;
    }
}

static void AVX2XorMemorySources(
    void* destination,
    const void* initial_source,
    const void* const* sources,
    uint32_t source_count,
    uint64_t byte_count)
{
    const void* waiting[8];
    unsigned waiting_count = 0;
    const void* accumulator = initial_source;
    bool wrote_destination = false;
    for (uint32_t i = 0; i < source_count; ++i)
    {
        if (!sources[i])
            continue;
        waiting[waiting_count++] = sources[i];
        if (waiting_count == 8)
        {
            AVX2XorMemorySourceGroup<8>(
                destination, accumulator, waiting, byte_count);
            accumulator = destination;
            wrote_destination = true;
            waiting_count = 0;
        }
    }
    if (waiting_count == 7)
    {
        AVX2XorMemorySourceGroup<7>(
            destination, accumulator, waiting, byte_count);
        return;
    }
    if (!wrote_destination)
        std::memcpy(destination, initial_source, static_cast<size_t>(byte_count));
    while (waiting_count >= 2)
    {
        AVX2XorMemory2To1(
            destination, waiting[0], waiting[1], byte_count);
        waiting_count -= 2;
        for (unsigned i = 0; i < waiting_count; ++i)
            waiting[i] = waiting[i + 2];
    }
    if (waiting_count != 0)
        AVX2XorMemory(destination, waiting[0], byte_count);
}

#else
static void AVX2XorMemorySources(
    void* destination,
    const void* initial_source,
    const void* const* sources,
    uint32_t source_count,
    uint64_t byte_count)
{
    std::memcpy(destination, initial_source, static_cast<size_t>(byte_count));
    const void* waiting[6];
    unsigned waiting_count = 0;
    for (uint32_t i = 0; i < source_count; ++i)
    {
        if (!sources[i])
            continue;
        waiting[waiting_count++] = sources[i];
        if (waiting_count == 6)
        {
            uint8_t* output = static_cast<uint8_t*>(destination);
            const uint8_t* input[6] = {
                static_cast<const uint8_t*>(waiting[0]),
                static_cast<const uint8_t*>(waiting[1]),
                static_cast<const uint8_t*>(waiting[2]),
                static_cast<const uint8_t*>(waiting[3]),
                static_cast<const uint8_t*>(waiting[4]),
                static_cast<const uint8_t*>(waiting[5])
            };
            uint64_t offset = 0;
            while (byte_count - offset >= 32)
            {
                __m256i result = _mm256_loadu_si256(
                    reinterpret_cast<const __m256i*>(output + offset));
                for (unsigned lane = 0; lane < 6; ++lane)
                    result = _mm256_xor_si256(result, _mm256_loadu_si256(
                        reinterpret_cast<const __m256i*>(
                            input[lane] + offset)));
                _mm256_storeu_si256(
                    reinterpret_cast<__m256i*>(output + offset), result);
                offset += 32;
            }
            while (offset < byte_count)
            {
                output[offset] ^= input[0][offset] ^ input[1][offset] ^
                    input[2][offset] ^ input[3][offset] ^ input[4][offset] ^
                    input[5][offset];
                ++offset;
            }
            waiting_count = 0;
        }
    }
    while (waiting_count >= 2)
    {
        AVX2XorMemory2To1(
            destination, waiting[0], waiting[1], byte_count);
        waiting_count -= 2;
        for (unsigned i = 0; i < waiting_count; ++i)
            waiting[i] = waiting[i + 2];
    }
    if (waiting_count != 0)
        AVX2XorMemory(destination, waiting[0], byte_count);
}
#endif

#if !defined(LEO2_AVX512_VARIANT)
static void AVX2XorMemoryDense(
    void* destination,
    const void* const* sources,
    uint32_t source_count,
    uint64_t byte_count)
{
    if (source_count == 2)
    {
        AVX2XorMemoryDenseGroup<2>(
            destination, sources, byte_count);
        return;
    }
    if (source_count == 4)
    {
        AVX2XorMemoryDenseGroup<4>(
            destination, sources, byte_count);
        return;
    }
    LEO_DEBUG_ASSERT(false);
}
#endif

static void AVX2XorMemory4(
    void* destination0, const void* source0,
    void* destination1, const void* source1,
    void* destination2, const void* source2,
    void* destination3, const void* source3,
    uint64_t byte_count)
{
    uint8_t* output0 = static_cast<uint8_t*>(destination0);
    uint8_t* output1 = static_cast<uint8_t*>(destination1);
    uint8_t* output2 = static_cast<uint8_t*>(destination2);
    uint8_t* output3 = static_cast<uint8_t*>(destination3);
    const uint8_t* input0 = static_cast<const uint8_t*>(source0);
    const uint8_t* input1 = static_cast<const uint8_t*>(source1);
    const uint8_t* input2 = static_cast<const uint8_t*>(source2);
    const uint8_t* input3 = static_cast<const uint8_t*>(source3);
    while (byte_count >= 32)
    {
        const __m256i result0 = _mm256_xor_si256(
            _mm256_loadu_si256(reinterpret_cast<const __m256i*>(output0)),
            _mm256_loadu_si256(reinterpret_cast<const __m256i*>(input0)));
        const __m256i result1 = _mm256_xor_si256(
            _mm256_loadu_si256(reinterpret_cast<const __m256i*>(output1)),
            _mm256_loadu_si256(reinterpret_cast<const __m256i*>(input1)));
        const __m256i result2 = _mm256_xor_si256(
            _mm256_loadu_si256(reinterpret_cast<const __m256i*>(output2)),
            _mm256_loadu_si256(reinterpret_cast<const __m256i*>(input2)));
        const __m256i result3 = _mm256_xor_si256(
            _mm256_loadu_si256(reinterpret_cast<const __m256i*>(output3)),
            _mm256_loadu_si256(reinterpret_cast<const __m256i*>(input3)));
        _mm256_storeu_si256(reinterpret_cast<__m256i*>(output0), result0);
        _mm256_storeu_si256(reinterpret_cast<__m256i*>(output1), result1);
        _mm256_storeu_si256(reinterpret_cast<__m256i*>(output2), result2);
        _mm256_storeu_si256(reinterpret_cast<__m256i*>(output3), result3);
        output0 += 32;
        output1 += 32;
        output2 += 32;
        output3 += 32;
        input0 += 32;
        input1 += 32;
        input2 += 32;
        input3 += 32;
        byte_count -= 32;
    }
    while (byte_count-- != 0)
    {
        *output0++ ^= *input0++;
        *output1++ ^= *input1++;
        *output2++ ^= *input2++;
        *output3++ ^= *input3++;
    }
}

#if defined(_MSC_VER)
#define LEO2_AVX2_FORCE_INLINE __forceinline
#else
#define LEO2_AVX2_FORCE_INLINE inline __attribute__((always_inline))
#endif

#ifdef LEO_HAS_FF8
static void AVX2FF8IFFTButterfly4Nonzero(
    void* value0_pointer, void* value1_pointer,
    void* value2_pointer, void* value3_pointer,
    uint16_t log01, uint16_t log23, uint16_t log02,
    uint64_t byte_count)
{
    uint8_t* value0 = static_cast<uint8_t*>(value0_pointer);
    uint8_t* value1 = static_cast<uint8_t*>(value1_pointer);
    uint8_t* value2 = static_cast<uint8_t*>(value2_pointer);
    uint8_t* value3 = static_cast<uint8_t*>(value3_pointer);
    const __m256i low01 = BroadcastTable(FF8Tables[log01].low);
    const __m256i high01 = BroadcastTable(FF8Tables[log01].high);
    const __m256i low23 = BroadcastTable(FF8Tables[log23].low);
    const __m256i high23 = BroadcastTable(FF8Tables[log23].high);
    const __m256i low02 = BroadcastTable(FF8Tables[log02].low);
    const __m256i high02 = BroadcastTable(FF8Tables[log02].high);

    while (byte_count >= 32)
    {
        __m256i x0 = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(value0));
        __m256i x1 = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(value1));
        __m256i x2 = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(value2));
        __m256i x3 = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(value3));
        x1 = _mm256_xor_si256(x1, x0);
        x0 = _mm256_xor_si256(x0,
            AVX2FF8ProductVector(x1, low01, high01));
        x3 = _mm256_xor_si256(x3, x2);
        x2 = _mm256_xor_si256(x2,
            AVX2FF8ProductVector(x3, low23, high23));
        x2 = _mm256_xor_si256(x2, x0);
        x3 = _mm256_xor_si256(x3, x1);
        x0 = _mm256_xor_si256(x0,
            AVX2FF8ProductVector(x2, low02, high02));
        x1 = _mm256_xor_si256(x1,
            AVX2FF8ProductVector(x3, low02, high02));
        _mm256_storeu_si256(reinterpret_cast<__m256i*>(value0), x0);
        _mm256_storeu_si256(reinterpret_cast<__m256i*>(value1), x1);
        _mm256_storeu_si256(reinterpret_cast<__m256i*>(value2), x2);
        _mm256_storeu_si256(reinterpret_cast<__m256i*>(value3), x3);
        value0 += 32;
        value1 += 32;
        value2 += 32;
        value3 += 32;
        byte_count -= 32;
    }
    while (byte_count-- != 0)
    {
        uint8_t x0 = *value0;
        uint8_t x1 = *value1;
        uint8_t x2 = *value2;
        uint8_t x3 = *value3;
        x1 ^= x0;
        x0 ^= FF8Product(log01, x1);
        x3 ^= x2;
        x2 ^= FF8Product(log23, x3);
        x2 ^= x0;
        x3 ^= x1;
        x0 ^= FF8Product(log02, x2);
        x1 ^= FF8Product(log02, x3);
        *value0++ = x0;
        *value1++ = x1;
        *value2++ = x2;
        *value3++ = x3;
    }
}

static void AVX2FF8IFFTButterfly4Kernel(
    void* value0_pointer, void* value1_pointer,
    void* value2_pointer, void* value3_pointer,
    uint16_t log01, uint16_t log23, uint16_t log02,
    uint64_t byte_count)
{
    static const uint16_t kZeroSkew = 255;
    if (log01 != kZeroSkew && log23 != kZeroSkew && log02 != kZeroSkew)
    {
        AVX2FF8IFFTButterfly4Nonzero(
            value0_pointer, value1_pointer, value2_pointer, value3_pointer,
            log01, log23, log02, byte_count);
        return;
    }

    uint8_t* value0 = static_cast<uint8_t*>(value0_pointer);
    uint8_t* value1 = static_cast<uint8_t*>(value1_pointer);
    uint8_t* value2 = static_cast<uint8_t*>(value2_pointer);
    uint8_t* value3 = static_cast<uint8_t*>(value3_pointer);

    __m256i low01 = _mm256_setzero_si256();
    __m256i high01 = _mm256_setzero_si256();
    __m256i low23 = _mm256_setzero_si256();
    __m256i high23 = _mm256_setzero_si256();
    __m256i low02 = _mm256_setzero_si256();
    __m256i high02 = _mm256_setzero_si256();
    if (log01 != kZeroSkew)
    {
        low01 = BroadcastTable(FF8Tables[log01].low);
        high01 = BroadcastTable(FF8Tables[log01].high);
    }
    if (log23 != kZeroSkew)
    {
        low23 = BroadcastTable(FF8Tables[log23].low);
        high23 = BroadcastTable(FF8Tables[log23].high);
    }
    if (log02 != kZeroSkew)
    {
        low02 = BroadcastTable(FF8Tables[log02].low);
        high02 = BroadcastTable(FF8Tables[log02].high);
    }

    while (byte_count >= 32)
    {
        __m256i x0 = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(value0));
        __m256i x1 = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(value1));
        __m256i x2 = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(value2));
        __m256i x3 = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(value3));
        x1 = _mm256_xor_si256(x1, x0);
        if (log01 != kZeroSkew)
            x0 = _mm256_xor_si256(x0,
                AVX2FF8ProductVector(x1, low01, high01));
        x3 = _mm256_xor_si256(x3, x2);
        if (log23 != kZeroSkew)
            x2 = _mm256_xor_si256(x2,
                AVX2FF8ProductVector(x3, low23, high23));
        x2 = _mm256_xor_si256(x2, x0);
        x3 = _mm256_xor_si256(x3, x1);
        if (log02 != kZeroSkew)
        {
            x0 = _mm256_xor_si256(x0,
                AVX2FF8ProductVector(x2, low02, high02));
            x1 = _mm256_xor_si256(x1,
                AVX2FF8ProductVector(x3, low02, high02));
        }
        _mm256_storeu_si256(reinterpret_cast<__m256i*>(value0), x0);
        _mm256_storeu_si256(reinterpret_cast<__m256i*>(value1), x1);
        _mm256_storeu_si256(reinterpret_cast<__m256i*>(value2), x2);
        _mm256_storeu_si256(reinterpret_cast<__m256i*>(value3), x3);
        value0 += 32;
        value1 += 32;
        value2 += 32;
        value3 += 32;
        byte_count -= 32;
    }
    while (byte_count-- != 0)
    {
        uint8_t x0 = *value0;
        uint8_t x1 = *value1;
        uint8_t x2 = *value2;
        uint8_t x3 = *value3;
        x1 ^= x0;
        if (log01 != kZeroSkew)
            x0 ^= FF8Product(log01, x1);
        x3 ^= x2;
        if (log23 != kZeroSkew)
            x2 ^= FF8Product(log23, x3);
        x2 ^= x0;
        x3 ^= x1;
        if (log02 != kZeroSkew)
        {
            x0 ^= FF8Product(log02, x2);
            x1 ^= FF8Product(log02, x3);
        }
        *value0++ = x0;
        *value1++ = x1;
        *value2++ = x2;
        *value3++ = x3;
    }
}

// The final radix-four layer for every message block after the first feeds an
// XOR accumulator and never observes the temporary work again.  Keep the
// transformed values in registers and accumulate them directly, matching the
// legacy encoder's single memory pass instead of storing work and rereading it
// through AVX2XorMemory4.
template<bool AllNonzero, bool Input3Zero = false>
static void AVX2FF8IFFTButterfly4XorKernel(
    const void* value0_pointer, const void* value1_pointer,
    const void* value2_pointer, const void* value3_pointer,
    void* output0_pointer, void* output1_pointer,
    void* output2_pointer, void* output3_pointer,
    uint16_t log01, uint16_t log23, uint16_t log02,
    uint64_t byte_count)
{
    static const uint16_t kZeroSkew = 255;
    const uint8_t* value0 = static_cast<const uint8_t*>(value0_pointer);
    const uint8_t* value1 = static_cast<const uint8_t*>(value1_pointer);
    const uint8_t* value2 = static_cast<const uint8_t*>(value2_pointer);
    const uint8_t* value3 = static_cast<const uint8_t*>(value3_pointer);
    uint8_t* output0 = static_cast<uint8_t*>(output0_pointer);
    uint8_t* output1 = static_cast<uint8_t*>(output1_pointer);
    uint8_t* output2 = static_cast<uint8_t*>(output2_pointer);
    uint8_t* output3 = static_cast<uint8_t*>(output3_pointer);

    __m256i low01 = _mm256_setzero_si256();
    __m256i high01 = _mm256_setzero_si256();
    __m256i low23 = _mm256_setzero_si256();
    __m256i high23 = _mm256_setzero_si256();
    __m256i low02 = _mm256_setzero_si256();
    __m256i high02 = _mm256_setzero_si256();
    if (AllNonzero || log01 != kZeroSkew)
    {
        low01 = BroadcastTable(FF8Tables[log01].low);
        high01 = BroadcastTable(FF8Tables[log01].high);
    }
    if (AllNonzero || log23 != kZeroSkew)
    {
        low23 = BroadcastTable(FF8Tables[log23].low);
        high23 = BroadcastTable(FF8Tables[log23].high);
    }
    if (AllNonzero || log02 != kZeroSkew)
    {
        low02 = BroadcastTable(FF8Tables[log02].low);
        high02 = BroadcastTable(FF8Tables[log02].high);
    }

    while (byte_count >= 32)
    {
        __m256i x0 = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(value0));
        __m256i x1 = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(value1));
        __m256i x2 = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(value2));
        __m256i x3 = Input3Zero
            ? _mm256_setzero_si256()
            : _mm256_loadu_si256(
                reinterpret_cast<const __m256i*>(value3));
        x1 = _mm256_xor_si256(x1, x0);
        if (AllNonzero || log01 != kZeroSkew)
            x0 = _mm256_xor_si256(x0,
                AVX2FF8ProductVector(x1, low01, high01));
        x3 = _mm256_xor_si256(x3, x2);
        if (AllNonzero || log23 != kZeroSkew)
            x2 = _mm256_xor_si256(x2,
                AVX2FF8ProductVector(x3, low23, high23));
        x2 = _mm256_xor_si256(x2, x0);
        x3 = _mm256_xor_si256(x3, x1);
        if (AllNonzero || log02 != kZeroSkew)
        {
            x0 = _mm256_xor_si256(x0,
                AVX2FF8ProductVector(x2, low02, high02));
            x1 = _mm256_xor_si256(x1,
                AVX2FF8ProductVector(x3, low02, high02));
        }
        x0 = _mm256_xor_si256(x0, _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(output0)));
        x1 = _mm256_xor_si256(x1, _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(output1)));
        x2 = _mm256_xor_si256(x2, _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(output2)));
        x3 = _mm256_xor_si256(x3, _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(output3)));
        _mm256_storeu_si256(reinterpret_cast<__m256i*>(output0), x0);
        _mm256_storeu_si256(reinterpret_cast<__m256i*>(output1), x1);
        _mm256_storeu_si256(reinterpret_cast<__m256i*>(output2), x2);
        _mm256_storeu_si256(reinterpret_cast<__m256i*>(output3), x3);
        value0 += 32;
        value1 += 32;
        value2 += 32;
        if (!Input3Zero)
            value3 += 32;
        output0 += 32;
        output1 += 32;
        output2 += 32;
        output3 += 32;
        byte_count -= 32;
    }
    while (byte_count-- != 0)
    {
        uint8_t x0 = *value0++;
        uint8_t x1 = *value1++;
        uint8_t x2 = *value2++;
        uint8_t x3 = Input3Zero ? 0 : *value3++;
        x1 ^= x0;
        if (AllNonzero || log01 != kZeroSkew)
            x0 ^= FF8Product(log01, x1);
        x3 ^= x2;
        if (AllNonzero || log23 != kZeroSkew)
            x2 ^= FF8Product(log23, x3);
        x2 ^= x0;
        x3 ^= x1;
        if (AllNonzero || log02 != kZeroSkew)
        {
            x0 ^= FF8Product(log02, x2);
            x1 ^= FF8Product(log02, x3);
        }
        *output0++ ^= x0;
        *output1++ ^= x1;
        *output2++ ^= x2;
        *output3++ ^= x3;
    }
}

static void AVX2FF8IFFTButterfly4(
    void* value0, void* value1, void* value2, void* value3,
    uint16_t log01, uint16_t log23, uint16_t log02,
    uint64_t byte_count)
{
    AVX2FF8IFFTButterfly4Kernel(value0, value1, value2, value3,
        log01, log23, log02, byte_count);
}

static void AVX2FF8FFTButterfly4Fused(
    void* value0_pointer, void* value1_pointer,
    void* value2_pointer, void* value3_pointer,
    uint16_t log01, uint16_t log23, uint16_t log02,
    uint64_t byte_count)
{
    static const uint16_t kZeroSkew = 255;
    uint8_t* value0 = static_cast<uint8_t*>(value0_pointer);
    uint8_t* value1 = static_cast<uint8_t*>(value1_pointer);
    uint8_t* value2 = static_cast<uint8_t*>(value2_pointer);
    uint8_t* value3 = static_cast<uint8_t*>(value3_pointer);
    __m256i low01 = _mm256_setzero_si256();
    __m256i high01 = _mm256_setzero_si256();
    __m256i low23 = _mm256_setzero_si256();
    __m256i high23 = _mm256_setzero_si256();
    __m256i low02 = _mm256_setzero_si256();
    __m256i high02 = _mm256_setzero_si256();
    if (log01 != kZeroSkew)
    {
        low01 = BroadcastTable(FF8Tables[log01].low);
        high01 = BroadcastTable(FF8Tables[log01].high);
    }
    if (log23 != kZeroSkew)
    {
        low23 = BroadcastTable(FF8Tables[log23].low);
        high23 = BroadcastTable(FF8Tables[log23].high);
    }
    if (log02 != kZeroSkew)
    {
        low02 = BroadcastTable(FF8Tables[log02].low);
        high02 = BroadcastTable(FF8Tables[log02].high);
    }

    while (byte_count >= 32)
    {
        __m256i x0 = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(value0));
        __m256i x1 = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(value1));
        __m256i x2 = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(value2));
        __m256i x3 = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(value3));
        if (log02 != kZeroSkew)
        {
            x0 = _mm256_xor_si256(x0,
                AVX2FF8ProductVector(x2, low02, high02));
            x1 = _mm256_xor_si256(x1,
                AVX2FF8ProductVector(x3, low02, high02));
        }
        x2 = _mm256_xor_si256(x2, x0);
        x3 = _mm256_xor_si256(x3, x1);
        if (log01 != kZeroSkew)
            x0 = _mm256_xor_si256(x0,
                AVX2FF8ProductVector(x1, low01, high01));
        x1 = _mm256_xor_si256(x1, x0);
        if (log23 != kZeroSkew)
            x2 = _mm256_xor_si256(x2,
                AVX2FF8ProductVector(x3, low23, high23));
        x3 = _mm256_xor_si256(x3, x2);
        _mm256_storeu_si256(reinterpret_cast<__m256i*>(value0), x0);
        _mm256_storeu_si256(reinterpret_cast<__m256i*>(value1), x1);
        _mm256_storeu_si256(reinterpret_cast<__m256i*>(value2), x2);
        _mm256_storeu_si256(reinterpret_cast<__m256i*>(value3), x3);
        value0 += 32;
        value1 += 32;
        value2 += 32;
        value3 += 32;
        byte_count -= 32;
    }
    while (byte_count-- != 0)
    {
        uint8_t x0 = *value0;
        uint8_t x1 = *value1;
        uint8_t x2 = *value2;
        uint8_t x3 = *value3;
        if (log02 != kZeroSkew)
        {
            x0 ^= FF8Product(log02, x2);
            x1 ^= FF8Product(log02, x3);
        }
        x2 ^= x0;
        x3 ^= x1;
        if (log01 != kZeroSkew)
            x0 ^= FF8Product(log01, x1);
        x1 ^= x0;
        if (log23 != kZeroSkew)
            x2 ^= FF8Product(log23, x3);
        x3 ^= x2;
        *value0++ = x0;
        *value1++ = x1;
        *value2++ = x2;
        *value3++ = x3;
    }
}

static void AVX2FF8FFTButterfly4(
    void* value0, void* value1, void* value2, void* value3,
    uint16_t log01, uint16_t log23, uint16_t log02,
    uint64_t byte_count)
{
    static const uint16_t kZeroSkew = 255;
    static const uint64_t kFusedByteLimit = 1024;
    if (byte_count <= kFusedByteLimit)
    {
        AVX2FF8FFTButterfly4Fused(
            value0, value1, value2, value3,
            log01, log23, log02, byte_count);
        return;
    }
    if (log02 == kZeroSkew)
    {
        AVX2XorMemory(value2, value0, byte_count);
        AVX2XorMemory(value3, value1, byte_count);
    }
    else
    {
        AVX2FF8Butterfly2<false>(value0, value2, log02, byte_count);
        AVX2FF8Butterfly2<false>(value1, value3, log02, byte_count);
    }
    if (log01 == kZeroSkew)
        AVX2XorMemory(value1, value0, byte_count);
    else
        AVX2FF8Butterfly2<false>(value0, value1, log01, byte_count);
    if (log23 == kZeroSkew)
        AVX2XorMemory(value3, value2, byte_count);
    else
        AVX2FF8Butterfly2<false>(value2, value3, log23, byte_count);
}

template<bool Inverse, uint8_t LiveMask = 15>
static void AVX2FF8Butterfly4Out(
    const void* input0_pointer, const void* input1_pointer,
    const void* input2_pointer, const void* input3_pointer,
    void* output0_pointer, void* output1_pointer,
    void* output2_pointer, void* output3_pointer,
    uint16_t log01, uint16_t log23, uint16_t log02,
    uint64_t byte_count)
{
    static const uint16_t kZeroSkew = 255;
    const uint8_t* input0 = static_cast<const uint8_t*>(input0_pointer);
    const uint8_t* input1 = static_cast<const uint8_t*>(input1_pointer);
    const uint8_t* input2 = static_cast<const uint8_t*>(input2_pointer);
    const uint8_t* input3 = static_cast<const uint8_t*>(input3_pointer);
    uint8_t* output0 = static_cast<uint8_t*>(output0_pointer);
    uint8_t* output1 = static_cast<uint8_t*>(output1_pointer);
    uint8_t* output2 = static_cast<uint8_t*>(output2_pointer);
    uint8_t* output3 = static_cast<uint8_t*>(output3_pointer);
    __m256i low01 = _mm256_setzero_si256();
    __m256i high01 = _mm256_setzero_si256();
    __m256i low23 = _mm256_setzero_si256();
    __m256i high23 = _mm256_setzero_si256();
    __m256i low02 = _mm256_setzero_si256();
    __m256i high02 = _mm256_setzero_si256();
    if (log01 != kZeroSkew)
    {
        low01 = BroadcastTable(FF8Tables[log01].low);
        high01 = BroadcastTable(FF8Tables[log01].high);
    }
    if (log23 != kZeroSkew)
    {
        low23 = BroadcastTable(FF8Tables[log23].low);
        high23 = BroadcastTable(FF8Tables[log23].high);
    }
    if (log02 != kZeroSkew)
    {
        low02 = BroadcastTable(FF8Tables[log02].low);
        high02 = BroadcastTable(FF8Tables[log02].high);
    }
    while (byte_count >= 32)
    {
        __m256i x0 = (LiveMask & 1U) != 0
            ? _mm256_loadu_si256(reinterpret_cast<const __m256i*>(input0))
            : _mm256_setzero_si256();
        __m256i x1 = (LiveMask & 2U) != 0
            ? _mm256_loadu_si256(reinterpret_cast<const __m256i*>(input1))
            : _mm256_setzero_si256();
        __m256i x2 = (LiveMask & 4U) != 0
            ? _mm256_loadu_si256(reinterpret_cast<const __m256i*>(input2))
            : _mm256_setzero_si256();
        __m256i x3 = (LiveMask & 8U) != 0
            ? _mm256_loadu_si256(reinterpret_cast<const __m256i*>(input3))
            : _mm256_setzero_si256();
        if (Inverse)
        {
            x1 = _mm256_xor_si256(x1, x0);
            if (log01 != kZeroSkew)
                x0 = _mm256_xor_si256(x0,
                    AVX2FF8ProductVector(x1, low01, high01));
            x3 = _mm256_xor_si256(x3, x2);
            if (log23 != kZeroSkew)
                x2 = _mm256_xor_si256(x2,
                    AVX2FF8ProductVector(x3, low23, high23));
            x2 = _mm256_xor_si256(x2, x0);
            x3 = _mm256_xor_si256(x3, x1);
            if (log02 != kZeroSkew)
            {
                x0 = _mm256_xor_si256(x0,
                    AVX2FF8ProductVector(x2, low02, high02));
                x1 = _mm256_xor_si256(x1,
                    AVX2FF8ProductVector(x3, low02, high02));
            }
        }
        else
        {
            if (log02 != kZeroSkew)
            {
                x0 = _mm256_xor_si256(x0,
                    AVX2FF8ProductVector(x2, low02, high02));
                x1 = _mm256_xor_si256(x1,
                    AVX2FF8ProductVector(x3, low02, high02));
            }
            x2 = _mm256_xor_si256(x2, x0);
            x3 = _mm256_xor_si256(x3, x1);
            if (log01 != kZeroSkew)
                x0 = _mm256_xor_si256(x0,
                    AVX2FF8ProductVector(x1, low01, high01));
            x1 = _mm256_xor_si256(x1, x0);
            if (log23 != kZeroSkew)
                x2 = _mm256_xor_si256(x2,
                    AVX2FF8ProductVector(x3, low23, high23));
            x3 = _mm256_xor_si256(x3, x2);
        }
        _mm256_storeu_si256(reinterpret_cast<__m256i*>(output0), x0);
        _mm256_storeu_si256(reinterpret_cast<__m256i*>(output1), x1);
        _mm256_storeu_si256(reinterpret_cast<__m256i*>(output2), x2);
        _mm256_storeu_si256(reinterpret_cast<__m256i*>(output3), x3);
        if ((LiveMask & 1U) != 0)
            input0 += 32;
        if ((LiveMask & 2U) != 0)
            input1 += 32;
        if ((LiveMask & 4U) != 0)
            input2 += 32;
        if ((LiveMask & 8U) != 0)
            input3 += 32;
        output0 += 32;
        output1 += 32;
        output2 += 32;
        output3 += 32;
        byte_count -= 32;
    }
    while (byte_count-- != 0)
    {
        uint8_t x0 = (LiveMask & 1U) != 0 ? *input0++ : 0;
        uint8_t x1 = (LiveMask & 2U) != 0 ? *input1++ : 0;
        uint8_t x2 = (LiveMask & 4U) != 0 ? *input2++ : 0;
        uint8_t x3 = (LiveMask & 8U) != 0 ? *input3++ : 0;
        if (Inverse)
        {
            x1 ^= x0;
            if (log01 != kZeroSkew)
                x0 ^= FF8Product(log01, x1);
            x3 ^= x2;
            if (log23 != kZeroSkew)
                x2 ^= FF8Product(log23, x3);
            x2 ^= x0;
            x3 ^= x1;
            if (log02 != kZeroSkew)
            {
                x0 ^= FF8Product(log02, x2);
                x1 ^= FF8Product(log02, x3);
            }
        }
        else
        {
            if (log02 != kZeroSkew)
            {
                x0 ^= FF8Product(log02, x2);
                x1 ^= FF8Product(log02, x3);
            }
            x2 ^= x0;
            x3 ^= x1;
            if (log01 != kZeroSkew)
                x0 ^= FF8Product(log01, x1);
            x1 ^= x0;
            if (log23 != kZeroSkew)
                x2 ^= FF8Product(log23, x3);
            x3 ^= x2;
        }
        *output0++ = x0;
        *output1++ = x1;
        *output2++ = x2;
        *output3++ = x3;
    }
}

static void AVX2FF8IFFTButterfly4Out(
    const void* input0, const void* input1,
    const void* input2, const void* input3,
    void* output0, void* output1, void* output2, void* output3,
    uint16_t log01, uint16_t log23, uint16_t log02,
    uint64_t byte_count)
{
    AVX2FF8Butterfly4Out<true>(
        input0, input1, input2, input3,
        output0, output1, output2, output3,
        log01, log23, log02, byte_count);
}

#if defined(LEO2_AVX512_VARIANT)

#if defined(_MSC_VER)
# define LEO2_AVX2_NOINLINE __declspec(noinline)
#elif defined(__GNUC__) || defined(__clang__)
# define LEO2_AVX2_NOINLINE __attribute__((noinline))
#else
# define LEO2_AVX2_NOINLINE
#endif

static LEO2_AVX2_NOINLINE void AVX2FF8IFFTButterfly4OutLastZero(
    const void* input0, const void* input1, const void* input2,
    void* output0, void* output1, void* output2, void* output3,
    uint16_t log01, uint16_t log23, uint16_t log02,
    uint64_t byte_count)
{
    AVX2FF8Butterfly4Out<true, 7>(
        input0, input1, input2, NULL,
        output0, output1, output2, output3,
        log01, log23, log02, byte_count);
}

#undef LEO2_AVX2_NOINLINE

#endif // LEO2_AVX512_VARIANT

static void AVX2FF8WeightedIFFTButterfly4(
    const void* input0_pointer, const void* input1_pointer,
    const void* input2_pointer, const void* input3_pointer,
    void* output0_pointer, void* output1_pointer,
    void* output2_pointer, void* output3_pointer,
    uint16_t weight_log0, uint16_t weight_log1,
    uint16_t weight_log2, uint16_t weight_log3,
    uint8_t live_mask,
    uint16_t log01, uint16_t log23, uint16_t log02,
    uint64_t byte_count)
{
    static const uint16_t kModulus = 255;
#if !defined(LEO2_AVX512_VARIANT)
    const bool identity_weights =
        (weight_log0 == 0 || weight_log0 == kModulus) &&
        (weight_log1 == 0 || weight_log1 == kModulus) &&
        (weight_log2 == 0 || weight_log2 == kModulus) &&
        (weight_log3 == 0 || weight_log3 == kModulus);
    if (identity_weights)
    {
        // Encoder ragged groups are live prefixes.  Dispatch once outside the
        // byte loop so their absent sources become compile-time zero vectors;
        // the general weighted/masked implementation remains authoritative
        // for every other locator pattern.
        switch (live_mask)
        {
        case 1:
            AVX2FF8Butterfly4Out<true, 1>(
                input0_pointer, NULL, NULL, NULL,
                output0_pointer, output1_pointer,
                output2_pointer, output3_pointer,
                log01, log23, log02, byte_count);
            return;
        case 3:
            AVX2FF8Butterfly4Out<true, 3>(
                input0_pointer, input1_pointer, NULL, NULL,
                output0_pointer, output1_pointer,
                output2_pointer, output3_pointer,
                log01, log23, log02, byte_count);
            return;
        case 7:
            AVX2FF8Butterfly4Out<true, 7>(
                input0_pointer, input1_pointer, input2_pointer, NULL,
                output0_pointer, output1_pointer,
                output2_pointer, output3_pointer,
                log01, log23, log02, byte_count);
            return;
        default:
            break;
        }
    }
#endif
    const uint8_t* inputs[4] = {
        static_cast<const uint8_t*>(input0_pointer),
        static_cast<const uint8_t*>(input1_pointer),
        static_cast<const uint8_t*>(input2_pointer),
        static_cast<const uint8_t*>(input3_pointer)
    };
    uint8_t* outputs[4] = {
        static_cast<uint8_t*>(output0_pointer),
        static_cast<uint8_t*>(output1_pointer),
        static_cast<uint8_t*>(output2_pointer),
        static_cast<uint8_t*>(output3_pointer)
    };
#ifdef LEO_DEBUG
    const void* alias_inputs[4] = {
        input0_pointer, input1_pointer, input2_pointer, input3_pointer
    };
    const void* alias_outputs[4] = {
        output0_pointer, output1_pointer, output2_pointer, output3_pointer
    };
    LEO_DEBUG_ASSERT(IsWeightedIFFTButterfly4AliasingValid(
        alias_inputs, alias_outputs, live_mask, byte_count));
#endif
    const uint16_t weights[4] = {
        weight_log0, weight_log1, weight_log2, weight_log3
    };
    __m256i low01 = _mm256_setzero_si256();
    __m256i high01 = _mm256_setzero_si256();
    __m256i low23 = _mm256_setzero_si256();
    __m256i high23 = _mm256_setzero_si256();
    __m256i low02 = _mm256_setzero_si256();
    __m256i high02 = _mm256_setzero_si256();
    if (log01 != kModulus)
    {
        low01 = BroadcastTable(FF8Tables[log01].low);
        high01 = BroadcastTable(FF8Tables[log01].high);
    }
    if (log23 != kModulus)
    {
        low23 = BroadcastTable(FF8Tables[log23].low);
        high23 = BroadcastTable(FF8Tables[log23].high);
    }
    if (log02 != kModulus)
    {
        low02 = BroadcastTable(FF8Tables[log02].low);
        high02 = BroadcastTable(FF8Tables[log02].high);
    }
    uint64_t offset = 0;
    while (byte_count - offset >= 32)
    {
        __m256i x[4];
        for (unsigned lane = 0; lane < 4; ++lane)
        {
            if ((live_mask & (1U << lane)) != 0)
            {
                x[lane] = _mm256_loadu_si256(
                    reinterpret_cast<const __m256i*>(inputs[lane] + offset));
                x[lane] = AVX2FF8ApplyWeight(x[lane], weights[lane]);
            }
            else
                x[lane] = _mm256_setzero_si256();
        }
        x[1] = _mm256_xor_si256(x[1], x[0]);
        if (log01 != kModulus)
            x[0] = _mm256_xor_si256(x[0],
                AVX2FF8ProductVector(x[1], low01, high01));
        x[3] = _mm256_xor_si256(x[3], x[2]);
        if (log23 != kModulus)
            x[2] = _mm256_xor_si256(x[2],
                AVX2FF8ProductVector(x[3], low23, high23));
        x[2] = _mm256_xor_si256(x[2], x[0]);
        x[3] = _mm256_xor_si256(x[3], x[1]);
        if (log02 != kModulus)
        {
            x[0] = _mm256_xor_si256(x[0],
                AVX2FF8ProductVector(x[2], low02, high02));
            x[1] = _mm256_xor_si256(x[1],
                AVX2FF8ProductVector(x[3], low02, high02));
        }
        for (unsigned lane = 0; lane < 4; ++lane)
            _mm256_storeu_si256(
                reinterpret_cast<__m256i*>(outputs[lane] + offset), x[lane]);
        offset += 32;
    }
    while (offset < byte_count)
    {
        uint8_t x[4] = {};
        for (unsigned lane = 0; lane < 4; ++lane)
        {
            if ((live_mask & (1U << lane)) == 0)
                continue;
            x[lane] = inputs[lane][offset];
            if (weights[lane] != 0 && weights[lane] != kModulus)
                x[lane] = FF8Product(weights[lane], x[lane]);
        }
        x[1] ^= x[0];
        if (log01 != kModulus)
            x[0] ^= FF8Product(log01, x[1]);
        x[3] ^= x[2];
        if (log23 != kModulus)
            x[2] ^= FF8Product(log23, x[3]);
        x[2] ^= x[0];
        x[3] ^= x[1];
        if (log02 != kModulus)
        {
            x[0] ^= FF8Product(log02, x[2]);
            x[1] ^= FF8Product(log02, x[3]);
        }
        for (unsigned lane = 0; lane < 4; ++lane)
            outputs[lane][offset] = x[lane];
        ++offset;
    }
}

static void AVX2FF8FFTButterfly4Out(
    const void* input0, const void* input1,
    const void* input2, const void* input3,
    void* output0, void* output1, void* output2, void* output3,
    uint16_t log01, uint16_t log23, uint16_t log02,
    uint64_t byte_count)
{
    AVX2FF8Butterfly4Out<false>(
        input0, input1, input2, input3,
        output0, output1, output2, output3,
        log01, log23, log02, byte_count);
}

template<bool Inverse>
static void AVX2FF8Butterfly4RangePrepared(
    void* const* work,
    uint32_t distance,
    uint16_t log01,
    uint16_t log23,
    uint16_t log02,
    uint64_t byte_count);

static void AVX2FF8IFFTButterfly4Range(
    void* const* work, unsigned distance,
    uint16_t log01, uint16_t log23, uint16_t log02,
    uint64_t byte_count, bool prefer_fused)
{
    (void)prefer_fused;
#if defined(LEO2_AVX512_VARIANT)
    for (unsigned i = 0; i < distance; ++i)
    {
        AVX2FF8IFFTButterfly4Kernel(
            work[i], work[i + distance],
            work[i + distance * 2U], work[i + distance * 3U],
            log01, log23, log02, byte_count);
    }
#else
    AVX2FF8Butterfly4RangePrepared<true>(
        work, distance, log01, log23, log02, byte_count);
#endif
}

static void AVX2FF8FFTButterfly4Range(
    void* const* work, unsigned distance,
    uint16_t log01, uint16_t log23, uint16_t log02,
    uint64_t byte_count, bool prefer_fused)
{
#if defined(LEO2_AVX512_VARIANT)
    if (prefer_fused)
    {
        for (unsigned i = 0; i < distance; ++i)
        {
            AVX2FF8FFTButterfly4Fused(
                work[i], work[i + distance],
                work[i + distance * 2U], work[i + distance * 3U],
                log01, log23, log02, byte_count);
        }
        return;
    }
#else
    static const uint64_t kFusedByteLimit = 1024;
    if (prefer_fused || byte_count <= kFusedByteLimit)
    {
        AVX2FF8Butterfly4RangePrepared<false>(
            work, distance, log01, log23, log02, byte_count);
        return;
    }
#endif
    // Retain the measured large-shard policy: each radix-two pass traverses
    // a complete shard before the next pass, rather than interleaving four
    // streams through the fused range kernel.
    for (unsigned i = 0; i < distance; ++i)
    {
        AVX2FF8FFTButterfly4(
            work[i], work[i + distance],
            work[i + distance * 2U], work[i + distance * 3U],
            log01, log23, log02, byte_count);
    }
}

static void AVX2FF8IFFTButterfly4XorRange(
    void* const* work, void* const* xor_output, unsigned distance,
    uint16_t log01, uint16_t log23, uint16_t log02,
    uint64_t byte_count)
{
    static const uint16_t kZeroSkew = 255;
    const bool all_nonzero =
        log01 != kZeroSkew && log23 != kZeroSkew && log02 != kZeroSkew;
    if (all_nonzero)
    {
        for (unsigned i = 0; i < distance; ++i)
            AVX2FF8IFFTButterfly4XorKernel<true>(
                work[i], work[i + distance],
                work[i + distance * 2U], work[i + distance * 3U],
                xor_output[i], xor_output[i + distance],
                xor_output[i + distance * 2U],
                xor_output[i + distance * 3U],
                log01, log23, log02, byte_count);
        return;
    }
    for (unsigned i = 0; i < distance; ++i)
        AVX2FF8IFFTButterfly4XorKernel<false>(
            work[i], work[i + distance],
            work[i + distance * 2U], work[i + distance * 3U],
            xor_output[i], xor_output[i + distance],
            xor_output[i + distance * 2U],
            xor_output[i + distance * 3U],
            log01, log23, log02, byte_count);
}

#if !defined(LEO2_AVX512_VARIANT)
struct AVX2FF8LinearTable
{
    __m256i low;
    __m256i high;
};

static LEO2_AVX2_FORCE_INLINE AVX2FF8LinearTable
AVX2FF8PrepareLinearTable(uint16_t first_log, uint16_t second_log)
{
    static const uint16_t kZeroSkew = 255;
    AVX2FF8LinearTable result;
    result.low = result.high = _mm256_setzero_si256();
    if (first_log != kZeroSkew)
    {
        result.low = BroadcastTable(FF8Tables[first_log].low);
        result.high = BroadcastTable(FF8Tables[first_log].high);
    }
    if (second_log != kZeroSkew)
    {
        result.low = _mm256_xor_si256(
            result.low, BroadcastTable(FF8Tables[second_log].low));
        result.high = _mm256_xor_si256(
            result.high, BroadcastTable(FF8Tables[second_log].high));
    }
    return result;
}

template<uint32_t OriginalCount>
static LEO2_AVX2_FORCE_INLINE void AVX2FF8HighEncodeT2Vector(
    const uint8_t* const* input,
    uint8_t* output0,
    uint8_t* output1,
    uint64_t offset,
    const AVX2FF8LinearTable& first,
    const AVX2FF8LinearTable& second)
{
    static_assert(OriginalCount == 2 || OriginalCount == 3,
        "T=2 direct encoder instantiated outside its generated K set");
    // Let a,b be the first message block, c the optional shortened-tail
    // source, u/v the two block-IFFT factors, and f the final FFT factor.
    // Expanding the two radix-two transforms over characteristic two gives:
    //
    //   q = (u + f) (a + b) + (v + f) c
    //   parity[0] = a + c + q, parity[1] = b + q.
    //
    // Addition of fixed field maps is XOR of their nibble tables, so this
    // executes the exact legacy transform in one pass with one product for
    // K=2 and two independent products for K=3.  The mature path otherwise
    // writes and rereads both accumulator shards for the final FFT, and K=3
    // also materializes a zero-padded second source block.
    const __m256i x0 = _mm256_loadu_si256(
        reinterpret_cast<const __m256i*>(input[0] + offset));
    const __m256i x1 = _mm256_loadu_si256(
        reinterpret_cast<const __m256i*>(input[1] + offset));
    const __m256i difference = _mm256_xor_si256(x0, x1);
    __m256i common = AVX2FF8ProductVector(
        difference, first.low, first.high);
    __m256i parity0 = x0;
    if (OriginalCount == 3)
    {
        const __m256i x2 = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(input[2] + offset));
        common = _mm256_xor_si256(common,
            AVX2FF8ProductVector(x2, second.low, second.high));
        parity0 = _mm256_xor_si256(parity0, x2);
    }
    parity0 = _mm256_xor_si256(parity0, common);
    const __m256i parity1 = _mm256_xor_si256(x1, common);
    _mm256_storeu_si256(
        reinterpret_cast<__m256i*>(output0 + offset), parity0);
    _mm256_storeu_si256(
        reinterpret_cast<__m256i*>(output1 + offset), parity1);
}

template<uint32_t OriginalCount>
#if defined(_MSC_VER)
# define LEO2_AVX2_T2_NOINLINE __declspec(noinline)
#elif defined(__GNUC__) || defined(__clang__)
# define LEO2_AVX2_T2_NOINLINE __attribute__((noinline))
#else
# define LEO2_AVX2_T2_NOINLINE
#endif
static LEO2_AVX2_T2_NOINLINE void AVX2FF8HighEncodeT2Direct(
    const void* const* data,
    void* const* work,
    const uint8_t* inverse_skew,
    const uint8_t* forward_skew,
    uint64_t byte_count)
{
    static const uint16_t kZeroSkew = 255;
    static_assert(OriginalCount == 2 || OriginalCount == 3,
        "T=2 direct encoder instantiated outside its generated K set");
    const AVX2FF8LinearTable first = AVX2FF8PrepareLinearTable(
        inverse_skew[1], forward_skew[1]);
    AVX2FF8LinearTable second;
    second.low = second.high = _mm256_setzero_si256();
    if (OriginalCount == 3)
    {
        second = AVX2FF8PrepareLinearTable(
            inverse_skew[3], forward_skew[1]);
    }

    const uint8_t* input[OriginalCount];
    for (uint32_t i = 0; i < OriginalCount; ++i)
        input[i] = static_cast<const uint8_t*>(data[i]);
    uint8_t* output0 = static_cast<uint8_t*>(work[0]);
    uint8_t* output1 = static_cast<uint8_t*>(work[1]);

    uint64_t offset = 0;
    for (; byte_count - offset >= 64U; offset += 64U)
    {
        AVX2FF8HighEncodeT2Vector<OriginalCount>(
            input, output0, output1, offset, first, second);
        AVX2FF8HighEncodeT2Vector<OriginalCount>(
            input, output0, output1, offset + 32U, first, second);
    }
    for (; byte_count - offset >= 32U; offset += 32U)
    {
        AVX2FF8HighEncodeT2Vector<OriginalCount>(
            input, output0, output1, offset, first, second);
    }
    for (; offset < byte_count; ++offset)
    {
        const uint8_t difference = static_cast<uint8_t>(
            input[0][offset] ^ input[1][offset]);
        uint8_t common = 0;
        if (inverse_skew[1] != kZeroSkew)
            common ^= FF8Product(inverse_skew[1], difference);
        if (forward_skew[1] != kZeroSkew)
            common ^= FF8Product(forward_skew[1], difference);
        uint8_t parity0 = input[0][offset];
        if (OriginalCount == 3)
        {
            const uint8_t x2 = input[2][offset];
            if (inverse_skew[3] != kZeroSkew)
                common ^= FF8Product(inverse_skew[3], x2);
            if (forward_skew[1] != kZeroSkew)
                common ^= FF8Product(forward_skew[1], x2);
            parity0 ^= x2;
        }
        output0[offset] = static_cast<uint8_t>(parity0 ^ common);
        output1[offset] = static_cast<uint8_t>(input[1][offset] ^ common);
    }
}
#undef LEO2_AVX2_T2_NOINLINE

struct AVX2FF8T4Tables
{
    __m256i low01;
    __m256i high01;
    __m256i low23;
    __m256i high23;
    __m256i low02;
    __m256i high02;
    bool has01;
    bool has23;
    bool has02;
};

static LEO2_AVX2_FORCE_INLINE AVX2FF8T4Tables AVX2FF8PrepareT4Tables(
    uint16_t log01, uint16_t log23, uint16_t log02)
{
    static const uint16_t kZeroSkew = 255;
    AVX2FF8T4Tables result;
    result.low01 = result.high01 = _mm256_setzero_si256();
    result.low23 = result.high23 = _mm256_setzero_si256();
    result.low02 = result.high02 = _mm256_setzero_si256();
    result.has01 = log01 != kZeroSkew;
    result.has23 = log23 != kZeroSkew;
    result.has02 = log02 != kZeroSkew;
    if (result.has01)
    {
        result.low01 = BroadcastTable(FF8Tables[log01].low);
        result.high01 = BroadcastTable(FF8Tables[log01].high);
    }
    if (result.has23)
    {
        result.low23 = BroadcastTable(FF8Tables[log23].low);
        result.high23 = BroadcastTable(FF8Tables[log23].high);
    }
    if (result.has02)
    {
        result.low02 = BroadcastTable(FF8Tables[log02].low);
        result.high02 = BroadcastTable(FF8Tables[log02].high);
    }
    return result;
}

static LEO2_AVX2_FORCE_INLINE void AVX2FF8T4Inverse(
    __m256i& x0, __m256i& x1, __m256i& x2, __m256i& x3,
    const AVX2FF8T4Tables& tables)
{
    x1 = _mm256_xor_si256(x1, x0);
    if (tables.has01)
        x0 = _mm256_xor_si256(x0,
            AVX2FF8ProductVector(x1, tables.low01, tables.high01));
    x3 = _mm256_xor_si256(x3, x2);
    if (tables.has23)
        x2 = _mm256_xor_si256(x2,
            AVX2FF8ProductVector(x3, tables.low23, tables.high23));
    x2 = _mm256_xor_si256(x2, x0);
    x3 = _mm256_xor_si256(x3, x1);
    if (tables.has02)
    {
        x0 = _mm256_xor_si256(x0,
            AVX2FF8ProductVector(x2, tables.low02, tables.high02));
        x1 = _mm256_xor_si256(x1,
            AVX2FF8ProductVector(x3, tables.low02, tables.high02));
    }
}

#if !defined(LEO2_GFNI_VARIANT)
static LEO2_AVX2_FORCE_INLINE __m256i AVX2LoadGF8Exact127Tail(
    const uint8_t* input)
{
    /*
        Assemble public bytes 96..126 plus the shortened zero byte without
        reading source[127].  The overlapping high load remains wholly inside
        the 127-byte object: it reads source[111..126], then drops byte 111.
    */
    const __m128i low = _mm_loadu_si128(
        reinterpret_cast<const __m128i*>(input + 96));
    const __m128i overlapping_high = _mm_loadu_si128(
        reinterpret_cast<const __m128i*>(input + 111));
    const __m128i high = _mm_srli_si128(overlapping_high, 1);
    return _mm256_inserti128_si256(
        _mm256_castsi128_si256(low), high, 1);
}

static LEO2_AVX2_FORCE_INLINE void AVX2FF8IFFTGroupExact127(
    const void* const* input_pointer,
    void* const* output_pointer,
    uint8_t live_mask,
    uint16_t multiplier_log01,
    uint16_t multiplier_log23,
    uint16_t multiplier_log02)
{
    const uint8_t* const input[4] = {
        static_cast<const uint8_t*>(input_pointer[0]),
        static_cast<const uint8_t*>(input_pointer[1]),
        static_cast<const uint8_t*>(input_pointer[2]),
        static_cast<const uint8_t*>(input_pointer[3])
    };
    uint8_t* const output[4] = {
        static_cast<uint8_t*>(output_pointer[0]),
        static_cast<uint8_t*>(output_pointer[1]),
        static_cast<uint8_t*>(output_pointer[2]),
        static_cast<uint8_t*>(output_pointer[3])
    };
    const AVX2FF8T4Tables tables = AVX2FF8PrepareT4Tables(
        multiplier_log01, multiplier_log23, multiplier_log02);

#if defined(__clang__)
# pragma clang loop unroll(disable)
#elif defined(__GNUC__)
# pragma GCC unroll 1
#endif
    for (unsigned offset = 0; offset < 128; offset += 32)
    {
        __m256i x0 = _mm256_setzero_si256();
        __m256i x1 = _mm256_setzero_si256();
        __m256i x2 = _mm256_setzero_si256();
        __m256i x3 = _mm256_setzero_si256();
        if (offset == 96)
        {
            if ((live_mask & 1U) != 0)
                x0 = AVX2LoadGF8Exact127Tail(input[0]);
            if ((live_mask & 2U) != 0)
                x1 = AVX2LoadGF8Exact127Tail(input[1]);
            if ((live_mask & 4U) != 0)
                x2 = AVX2LoadGF8Exact127Tail(input[2]);
            if ((live_mask & 8U) != 0)
                x3 = AVX2LoadGF8Exact127Tail(input[3]);
        }
        else
        {
            if ((live_mask & 1U) != 0)
                x0 = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(
                    input[0] + offset));
            if ((live_mask & 2U) != 0)
                x1 = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(
                    input[1] + offset));
            if ((live_mask & 4U) != 0)
                x2 = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(
                    input[2] + offset));
            if ((live_mask & 8U) != 0)
                x3 = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(
                    input[3] + offset));
        }
        AVX2FF8T4Inverse(x0, x1, x2, x3, tables);
        _mm256_storeu_si256(
            reinterpret_cast<__m256i*>(output[0] + offset), x0);
        _mm256_storeu_si256(
            reinterpret_cast<__m256i*>(output[1] + offset), x1);
        _mm256_storeu_si256(
            reinterpret_cast<__m256i*>(output[2] + offset), x2);
        _mm256_storeu_si256(
            reinterpret_cast<__m256i*>(output[3] + offset), x3);
    }
}

void AVX2FF8IFFTBlockFromSourcesExact127(
    const void* const* sources,
    void** work,
    uint32_t count,
    const uint8_t* skew_lut)
{
    LEO_DEBUG_ASSERT(count > 4 && (count & 3U) == 0);
    for (uint32_t r = 0; r < count; r += 4)
    {
        uint8_t live_mask = 0;
        for (unsigned lane = 0; lane < 4; ++lane)
            if (sources[r + lane])
                live_mask |= static_cast<uint8_t>(1U << lane);
        if (live_mask == 0)
        {
            for (unsigned lane = 0; lane < 4; ++lane)
                memset(work[r + lane], 0, 128);
            continue;
        }
        AVX2FF8IFFTGroupExact127(
            sources + r, work + r, live_mask,
            skew_lut[r + 1], skew_lut[r + 3], skew_lut[r + 2]);
    }
}
#endif

static LEO2_AVX2_FORCE_INLINE void AVX2FF8T4Forward(
    __m256i& x0, __m256i& x1, __m256i& x2, __m256i& x3,
    const AVX2FF8T4Tables& tables)
{
    if (tables.has02)
    {
        x0 = _mm256_xor_si256(x0,
            AVX2FF8ProductVector(x2, tables.low02, tables.high02));
        x1 = _mm256_xor_si256(x1,
            AVX2FF8ProductVector(x3, tables.low02, tables.high02));
    }
    x2 = _mm256_xor_si256(x2, x0);
    x3 = _mm256_xor_si256(x3, x1);
    if (tables.has01)
        x0 = _mm256_xor_si256(x0,
            AVX2FF8ProductVector(x1, tables.low01, tables.high01));
    x1 = _mm256_xor_si256(x1, x0);
    if (tables.has23)
        x2 = _mm256_xor_si256(x2,
            AVX2FF8ProductVector(x3, tables.low23, tables.high23));
    x3 = _mm256_xor_si256(x3, x2);
}

template<uint32_t OriginalCount, uint32_t RecoveryCount>
static LEO2_AVX2_FORCE_INLINE void AVX2FF8HighEncodeT4BlocksPrepared(
    const void* const* data,
    void* const* work,
    uint64_t byte_count,
    const AVX2FF8T4Tables& inverse0,
    const AVX2FF8T4Tables& inverse1,
    const AVX2FF8T4Tables& inverse2,
    const AVX2FF8T4Tables& inverse3,
    const AVX2FF8T4Tables& forward)
{
    static_assert(
        (OriginalCount >= 3 && OriginalCount <= 11) ||
        (OriginalCount >= 13 && OriginalCount <= 14),
        "T=4 fused encoder instantiated outside its supported K set");
    static_assert(
        RecoveryCount == 3 || RecoveryCount == 4,
        "T=4 batch encoder instantiated outside its parity prefix");
    LEO_DEBUG_ASSERT((byte_count & 31U) == 0);

    const uint8_t* input[OriginalCount];
    for (uint32_t i = 0; i < OriginalCount; ++i)
        input[i] = static_cast<const uint8_t*>(data[i]);
    uint8_t* output0 = static_cast<uint8_t*>(work[0]);
    uint8_t* output1 = static_cast<uint8_t*>(work[1]);
    uint8_t* output2 = static_cast<uint8_t*>(work[2]);
    uint8_t* output3 = RecoveryCount == 4
        ? static_cast<uint8_t*>(work[3]) : NULL;

    // Keep the transformed accumulator live across all source blocks and the
    // final forward transform.  This removes both the temporary ragged block
    // and every full-shard accumulator round trip from the selected tiny-T
    // shapes.  OriginalCount is a template argument so missing suffix lanes
    // fold to a zero register without a branch in the byte loop.
    for (uint64_t offset = 0; offset < byte_count; offset += 32)
    {
        __m256i x0 = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(input[0] + offset));
        __m256i x1 = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(input[1] + offset));
        __m256i x2 = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(input[2] + offset));
        __m256i x3 = OriginalCount >= 4
            ? _mm256_loadu_si256(
                reinterpret_cast<const __m256i*>(input[3] + offset))
            : _mm256_setzero_si256();
        AVX2FF8T4Inverse(x0, x1, x2, x3, inverse0);

        if (OriginalCount > 4)
        {
            __m256i y0 = _mm256_loadu_si256(
                reinterpret_cast<const __m256i*>(input[4] + offset));
            __m256i y1 = OriginalCount >= 6
                ? _mm256_loadu_si256(
                    reinterpret_cast<const __m256i*>(input[5] + offset))
                : _mm256_setzero_si256();
            __m256i y2 = OriginalCount >= 7
                ? _mm256_loadu_si256(
                    reinterpret_cast<const __m256i*>(input[6] + offset))
                : _mm256_setzero_si256();
            __m256i y3 = OriginalCount >= 8
                ? _mm256_loadu_si256(
                    reinterpret_cast<const __m256i*>(input[7] + offset))
                : _mm256_setzero_si256();
            AVX2FF8T4Inverse(y0, y1, y2, y3, inverse1);
            x0 = _mm256_xor_si256(x0, y0);
            x1 = _mm256_xor_si256(x1, y1);
            x2 = _mm256_xor_si256(x2, y2);
            x3 = _mm256_xor_si256(x3, y3);
        }

        if (OriginalCount > 8)
        {
            __m256i z0 = _mm256_loadu_si256(
                reinterpret_cast<const __m256i*>(input[8] + offset));
            __m256i z1 = OriginalCount >= 10
                ? _mm256_loadu_si256(
                    reinterpret_cast<const __m256i*>(input[9] + offset))
                : _mm256_setzero_si256();
            __m256i z2 = OriginalCount >= 11
                ? _mm256_loadu_si256(
                    reinterpret_cast<const __m256i*>(input[10] + offset))
                : _mm256_setzero_si256();
            __m256i z3 = OriginalCount >= 12
                ? _mm256_loadu_si256(
                    reinterpret_cast<const __m256i*>(input[11] + offset))
                : _mm256_setzero_si256();
            AVX2FF8T4Inverse(z0, z1, z2, z3, inverse2);
            x0 = _mm256_xor_si256(x0, z0);
            x1 = _mm256_xor_si256(x1, z1);
            x2 = _mm256_xor_si256(x2, z2);
            x3 = _mm256_xor_si256(x3, z3);
        }

        if (OriginalCount > 12)
        {
            __m256i w0 = _mm256_loadu_si256(
                reinterpret_cast<const __m256i*>(input[12] + offset));
            __m256i w1 = OriginalCount >= 14
                ? _mm256_loadu_si256(
                    reinterpret_cast<const __m256i*>(input[13] + offset))
                : _mm256_setzero_si256();
            __m256i w2 = OriginalCount >= 15
                ? _mm256_loadu_si256(
                    reinterpret_cast<const __m256i*>(input[14] + offset))
                : _mm256_setzero_si256();
            __m256i w3 = OriginalCount >= 16
                ? _mm256_loadu_si256(
                    reinterpret_cast<const __m256i*>(input[15] + offset))
                : _mm256_setzero_si256();
            AVX2FF8T4Inverse(w0, w1, w2, w3, inverse3);
            x0 = _mm256_xor_si256(x0, w0);
            x1 = _mm256_xor_si256(x1, w1);
            x2 = _mm256_xor_si256(x2, w2);
            x3 = _mm256_xor_si256(x3, w3);
        }

        AVX2FF8T4Forward(x0, x1, x2, x3, forward);
        _mm256_storeu_si256(
            reinterpret_cast<__m256i*>(output0 + offset), x0);
        _mm256_storeu_si256(
            reinterpret_cast<__m256i*>(output1 + offset), x1);
        _mm256_storeu_si256(
            reinterpret_cast<__m256i*>(output2 + offset), x2);
        if (RecoveryCount == 4)
        {
            _mm256_storeu_si256(
                reinterpret_cast<__m256i*>(output3 + offset), x3);
        }
    }
}

template<uint32_t OriginalCount>
static void AVX2FF8HighEncodeT4Blocks(
    const void* const* data,
    void* const* work,
    const uint8_t* inverse_skew,
    const uint8_t* forward_skew,
    uint64_t byte_count)
{
    const AVX2FF8T4Tables inverse0 = AVX2FF8PrepareT4Tables(
        inverse_skew[1], inverse_skew[3], inverse_skew[2]);
    const AVX2FF8T4Tables forward = AVX2FF8PrepareT4Tables(
        forward_skew[1], forward_skew[3], forward_skew[2]);
    AVX2FF8T4Tables inverse1 = {};
    if (OriginalCount > 4)
        inverse1 = AVX2FF8PrepareT4Tables(
            inverse_skew[5], inverse_skew[7], inverse_skew[6]);
    AVX2FF8T4Tables inverse2 = {};
    if (OriginalCount > 8)
        inverse2 = AVX2FF8PrepareT4Tables(
            inverse_skew[9], inverse_skew[11], inverse_skew[10]);
    AVX2FF8T4Tables inverse3 = {};
    if (OriginalCount > 12)
        inverse3 = AVX2FF8PrepareT4Tables(
            inverse_skew[13], inverse_skew[15], inverse_skew[14]);
    AVX2FF8HighEncodeT4BlocksPrepared<OriginalCount, 4>(
        data, work, byte_count,
        inverse0, inverse1, inverse2, inverse3, forward);
}

#if !defined(LEO2_AVX512_VARIANT) && !defined(LEO2_GFNI_VARIANT)
template<uint32_t OriginalCount, uint32_t RecoveryCount>
static void AVX2FF8HighEncodeT4BatchBlocks(
    const void* const* data,
    void* const* recovery,
    uint32_t item_count,
    const uint8_t* inverse_skew,
    const uint8_t* forward_skew,
    uint64_t byte_count)
{
    const AVX2FF8T4Tables inverse0 = AVX2FF8PrepareT4Tables(
        inverse_skew[1], inverse_skew[3], inverse_skew[2]);
    const AVX2FF8T4Tables forward = AVX2FF8PrepareT4Tables(
        forward_skew[1], forward_skew[3], forward_skew[2]);
    AVX2FF8T4Tables inverse1 = {};
    if (OriginalCount > 4)
        inverse1 = AVX2FF8PrepareT4Tables(
            inverse_skew[5], inverse_skew[7], inverse_skew[6]);
    AVX2FF8T4Tables inverse2 = {};
    if (OriginalCount > 8)
        inverse2 = AVX2FF8PrepareT4Tables(
            inverse_skew[9], inverse_skew[11], inverse_skew[10]);
    AVX2FF8T4Tables inverse3 = {};
    if (OriginalCount > 12)
        inverse3 = AVX2FF8PrepareT4Tables(
            inverse_skew[13], inverse_skew[15], inverse_skew[14]);

    for (uint32_t item = 0; item < item_count; ++item)
    {
        AVX2FF8HighEncodeT4BlocksPrepared<
            OriginalCount, RecoveryCount>(
            data + static_cast<size_t>(item) * OriginalCount,
            recovery + static_cast<size_t>(item) * RecoveryCount,
            byte_count, inverse0, inverse1, inverse2, inverse3, forward);
    }
}

static void AVX2FF8HighEncodeT4Batch(
    const void* const* data,
    void* const* recovery,
    uint32_t item_count,
    uint32_t original_count,
    uint32_t recovery_count,
    const uint8_t* inverse_skew,
    const uint8_t* forward_skew,
    uint64_t byte_count)
{
    LEO_DEBUG_ASSERT(data != NULL);
    LEO_DEBUG_ASSERT(recovery != NULL);
    LEO_DEBUG_ASSERT(item_count != 0);
    LEO_DEBUG_ASSERT(recovery_count == 3 || recovery_count == 4);
    LEO_DEBUG_ASSERT(
        byte_count >= 32 && byte_count <= 16U * 1024U &&
        (byte_count & 31U) == 0);

    /* One exact tiny stripe does not amortize the batch wrapper.  Reuse the
       established register-resident circuit directly; multi-stripe bindings
       retain the prepared-table batch specializations below. */
    if (original_count == 4 && recovery_count == 4 &&
        item_count == 1 && byte_count == 64)
    {
        AVX2FF8HighEncodeT4Blocks<4>(
            data, recovery, inverse_skew, forward_skew, byte_count);
        return;
    }

    /* This callback is private to already-qualified Leopard2 codecs and its
       startup known-answer test.  Collapse the old recovery branch followed
       by an original-count jump table into one fixed-shape dispatch. */
#define LEO2_AVX2_T4_SHAPE_CASE(OriginalCount, RecoveryCount) \
    case (RecoveryCount << 4) | OriginalCount: \
        AVX2FF8HighEncodeT4BatchBlocks<OriginalCount, RecoveryCount>( \
            data, recovery, item_count, inverse_skew, forward_skew, \
            byte_count); \
        return

    switch ((recovery_count << 4) | original_count)
    {
    /* Retain the mature R=4-first instantiation order. */
    LEO2_AVX2_T4_SHAPE_CASE(3, 4);
    LEO2_AVX2_T4_SHAPE_CASE(4, 4);
    LEO2_AVX2_T4_SHAPE_CASE(5, 4);
    LEO2_AVX2_T4_SHAPE_CASE(6, 4);
    LEO2_AVX2_T4_SHAPE_CASE(7, 4);
    /* K=8 is selected only by the exact-byte packed public terminal.  The
       general ff8_high_encode_small route below deliberately still omits it. */
    LEO2_AVX2_T4_SHAPE_CASE(8, 4);
    LEO2_AVX2_T4_SHAPE_CASE(9, 4);
    LEO2_AVX2_T4_SHAPE_CASE(10, 4);
    LEO2_AVX2_T4_SHAPE_CASE(11, 4);
    LEO2_AVX2_T4_SHAPE_CASE(13, 4);
    LEO2_AVX2_T4_SHAPE_CASE(14, 4);
    LEO2_AVX2_T4_SHAPE_CASE(3, 3);
    LEO2_AVX2_T4_SHAPE_CASE(4, 3);
    LEO2_AVX2_T4_SHAPE_CASE(5, 3);
    LEO2_AVX2_T4_SHAPE_CASE(6, 3);
    LEO2_AVX2_T4_SHAPE_CASE(7, 3);
    LEO2_AVX2_T4_SHAPE_CASE(8, 3);
    LEO2_AVX2_T4_SHAPE_CASE(9, 3);
    LEO2_AVX2_T4_SHAPE_CASE(10, 3);
    LEO2_AVX2_T4_SHAPE_CASE(11, 3);
    LEO2_AVX2_T4_SHAPE_CASE(13, 3);
    LEO2_AVX2_T4_SHAPE_CASE(14, 3);
    default:
        LEO_DEBUG_BREAK;
        return;
    }

#undef LEO2_AVX2_T4_SHAPE_CASE
}

#endif

static void AVX2FF8HighEncodeSmall(
    const void* const* data,
    uint32_t original_count,
    void* const* work,
    uint32_t side,
    const uint8_t* inverse_skew,
    const uint8_t* forward_skew,
    uint64_t byte_count)
{
    static const uint16_t kZeroSkew = 255;

    // Below the established whole-transform thresholds, retaining the mature
    // inverse accumulation but fusing only the final forward T=4 butterfly
    // removes two extra full-shard memory passes without the register pressure
    // of the combined multi-block templates.
    const bool existing_fused_t4_threshold =
        // The backend qualification harness calls the callback directly with
        // sub-2-KiB buffers.  Production dispatch begins at 2 KiB, so retain
        // the register-fused implementation below that boundary to keep every
        // generated specialization covered by the startup known-answer test.
        byte_count < 2U * 1024U ||
        ((original_count == 5 || original_count == 6 ||
            original_count >= 9) && byte_count >= 2U * 1024U) ||
        ((original_count == 3 || original_count == 7) &&
            byte_count >= 4U * 1024U) ||
        (original_count == 4 && byte_count >= 8U * 1024U);
    const bool fused_t4_count = existing_fused_t4_threshold && (
        (original_count >= 3 && original_count <= 7) ||
        (original_count >= 9 && original_count <= 11)
#if !defined(LEO2_AVX512_VARIANT) && !defined(LEO2_GFNI_VARIANT)
        /* These four-block shapes are qualified only for the AVX2 backend.
           Avoid cloning their large fixed circuits into unrelated runtime
           backends until those variants have independent measurements. */
        || (original_count >= 13 && original_count <= 14)
#endif
        );
    if (side == 4 && fused_t4_count && (byte_count & 31U) == 0)
    {
        switch (original_count)
        {
        case 3: AVX2FF8HighEncodeT4Blocks<3>(
            data, work, inverse_skew, forward_skew, byte_count); return;
        case 4: AVX2FF8HighEncodeT4Blocks<4>(
            data, work, inverse_skew, forward_skew, byte_count); return;
        case 5: AVX2FF8HighEncodeT4Blocks<5>(
            data, work, inverse_skew, forward_skew, byte_count); return;
        case 6: AVX2FF8HighEncodeT4Blocks<6>(
            data, work, inverse_skew, forward_skew, byte_count); return;
        case 7: AVX2FF8HighEncodeT4Blocks<7>(
            data, work, inverse_skew, forward_skew, byte_count); return;
        case 9: AVX2FF8HighEncodeT4Blocks<9>(
            data, work, inverse_skew, forward_skew, byte_count); return;
        case 10: AVX2FF8HighEncodeT4Blocks<10>(
            data, work, inverse_skew, forward_skew, byte_count); return;
        case 11: AVX2FF8HighEncodeT4Blocks<11>(
            data, work, inverse_skew, forward_skew, byte_count); return;
#if !defined(LEO2_AVX512_VARIANT) && !defined(LEO2_GFNI_VARIANT)
        /* Reuse the already-emitted R=4 batch specializations instead of
           cloning another ~9 KiB pair for detached/sparse fallback layouts. */
        case 13: AVX2FF8HighEncodeT4BatchBlocks<13, 4>(
            data, work, 1, inverse_skew, forward_skew, byte_count); return;
        case 14: AVX2FF8HighEncodeT4BatchBlocks<14, 4>(
            data, work, 1, inverse_skew, forward_skew, byte_count); return;
#endif
        }
    }

    // Materialize the first complete message block directly into the
    // accumulator.  This avoids both a separate output clear and the
    // redundant zero-valued accumulator load in the first XOR transform.
    if (side == 2)
    {
        if (original_count == 2)
        {
            AVX2FF8HighEncodeT2Direct<2>(
                data, work, inverse_skew, forward_skew, byte_count);
            return;
        }
        if (original_count == 3)
        {
            AVX2FF8HighEncodeT2Direct<3>(
                data, work, inverse_skew, forward_skew, byte_count);
            return;
        }
        AVX2FF8IFFTButterfly2Out(
            data[0], data[1], work[0], work[1],
            inverse_skew[1], byte_count);
    }
    else if (original_count >= 4)
    {
        AVX2FF8IFFTButterfly4Out(
            data[0], data[1], data[2], data[3],
            work[0], work[1], work[2], work[3],
            inverse_skew[1], inverse_skew[3], inverse_skew[2],
            byte_count);
    }
    else
    {
        LEO_DEBUG_ASSERT(side == 4 && original_count == 3);
        AVX2FF8Butterfly4Out<true, 7>(
            data[0], data[1], data[2], NULL,
            work[0], work[1], work[2], work[3],
            inverse_skew[1], inverse_skew[3], inverse_skew[2],
            byte_count);
    }

    for (uint32_t base = side; base < original_count; base += side)
    {
        const uint32_t remaining = std::min(
            side, original_count - base);
        const void* source[4] = { NULL, NULL, NULL, NULL };
        const bool direct_last_zero = side == 4 && remaining == 3;
        if (remaining == side || direct_last_zero)
        {
            for (uint32_t lane = 0; lane < remaining; ++lane)
                source[lane] = data[base + lane];
        }
        else
        {
            // work[side..2*side) is the established encoder temporary half.
            // A high-rate block can have at most one shortened suffix, so the
            // only retained copy/zero pass is bounded independently of K.
            for (uint32_t lane = 0; lane < side; ++lane)
            {
                if (lane < remaining)
                    std::memcpy(work[side + lane], data[base + lane],
                        static_cast<size_t>(byte_count));
                else
                    std::memset(work[side + lane], 0,
                        static_cast<size_t>(byte_count));
                source[lane] = work[side + lane];
            }
        }

        if (side == 2)
        {
            const uint16_t log = inverse_skew[base + 1U];
            if (log == kZeroSkew)
            {
                AVX2XorMemory(work[0], source[0], byte_count);
                AVX2XorMemory2To1(
                    work[1], source[0], source[1], byte_count);
            }
            else
            {
                AVX2FF8IFFTButterfly2Xor(
                    source[0], source[1], work[0], work[1],
                    log, byte_count);
            }
        }
        else
        {
            const uint16_t log01 = inverse_skew[base + 1U];
            const uint16_t log23 = inverse_skew[base + 3U];
            const uint16_t log02 = inverse_skew[base + 2U];
            const bool all_nonzero = log01 != kZeroSkew &&
                log23 != kZeroSkew && log02 != kZeroSkew;
            if (direct_last_zero)
            {
                if (all_nonzero)
                    AVX2FF8IFFTButterfly4XorKernel<true, true>(
                        source[0], source[1], source[2], NULL,
                        work[0], work[1], work[2], work[3],
                        log01, log23, log02, byte_count);
                else
                    AVX2FF8IFFTButterfly4XorKernel<false, true>(
                        source[0], source[1], source[2], NULL,
                        work[0], work[1], work[2], work[3],
                        log01, log23, log02, byte_count);
            }
            else if (all_nonzero)
            {
                AVX2FF8IFFTButterfly4XorKernel<true>(
                    source[0], source[1], source[2], source[3],
                    work[0], work[1], work[2], work[3],
                    log01, log23, log02, byte_count);
            }
            else
                AVX2FF8IFFTButterfly4XorKernel<false>(
                    source[0], source[1], source[2], source[3],
                    work[0], work[1], work[2], work[3],
                    log01, log23, log02, byte_count);
        }
    }

    if (side == 2)
    {
        const uint16_t log = forward_skew[1];
        if (log == kZeroSkew)
            AVX2XorMemory(work[1], work[0], byte_count);
        else
            AVX2FF8FFTButterfly2(work[0], work[1], log, byte_count);
    }
    else
    {
        const bool final_only_fused_t4 =
            ((original_count == 3 || original_count == 7 ||
                original_count == 8) && byte_count < 4U * 1024U) ||
            (original_count == 4 && byte_count < 8U * 1024U);
        if (final_only_fused_t4)
        {
            AVX2FF8FFTButterfly4Fused(
                work[0], work[1], work[2], work[3],
                forward_skew[1], forward_skew[3], forward_skew[2],
                byte_count);
        }
        else
        {
            AVX2FF8FFTButterfly4(
                work[0], work[1], work[2], work[3],
                forward_skew[1], forward_skew[3], forward_skew[2],
                byte_count);
        }
    }
}
#endif

#if defined(LEO2_AVX512_VARIANT)

struct AVX2FF8T8Pair
{
    __m256i low;
    __m256i high;
};

static LEO2_AVX2_FORCE_INLINE void AVX2FF8T8Xor(
    AVX2FF8T8Pair& destination,
    const AVX2FF8T8Pair& source)
{
    destination.low = _mm256_xor_si256(destination.low, source.low);
    destination.high = _mm256_xor_si256(destination.high, source.high);
}

static LEO2_AVX2_FORCE_INLINE void AVX2FF8T8MulAdd(
    AVX2FF8T8Pair& destination,
    const AVX2FF8T8Pair& source,
    uint16_t log)
{
    const FF8NibbleTable& table = FF8Tables[log];
    const __m256i low_table = BroadcastTable(table.low);
    const __m256i high_table = BroadcastTable(table.high);
    destination.low = _mm256_xor_si256(destination.low,
        AVX2FF8ProductVector(source.low, low_table, high_table));
    destination.high = _mm256_xor_si256(destination.high,
        AVX2FF8ProductVector(source.high, low_table, high_table));
}

static LEO2_AVX2_FORCE_INLINE void AVX2FF8T8MulAdd2(
    AVX2FF8T8Pair& destination0,
    const AVX2FF8T8Pair& source0,
    AVX2FF8T8Pair& destination1,
    const AVX2FF8T8Pair& source1,
    uint16_t log)
{
    const FF8NibbleTable& table = FF8Tables[log];
    const __m256i low_table = BroadcastTable(table.low);
    const __m256i high_table = BroadcastTable(table.high);
    destination0.low = _mm256_xor_si256(destination0.low,
        AVX2FF8ProductVector(source0.low, low_table, high_table));
    destination0.high = _mm256_xor_si256(destination0.high,
        AVX2FF8ProductVector(source0.high, low_table, high_table));
    destination1.low = _mm256_xor_si256(destination1.low,
        AVX2FF8ProductVector(source1.low, low_table, high_table));
    destination1.high = _mm256_xor_si256(destination1.high,
        AVX2FF8ProductVector(source1.high, low_table, high_table));
}

static LEO2_AVX2_FORCE_INLINE void AVX2FF8T8MulAdd4(
    AVX2FF8T8Pair& destination0,
    const AVX2FF8T8Pair& source0,
    AVX2FF8T8Pair& destination1,
    const AVX2FF8T8Pair& source1,
    AVX2FF8T8Pair& destination2,
    const AVX2FF8T8Pair& source2,
    AVX2FF8T8Pair& destination3,
    const AVX2FF8T8Pair& source3,
    uint16_t log)
{
    const FF8NibbleTable& table = FF8Tables[log];
    const __m256i low_table = BroadcastTable(table.low);
    const __m256i high_table = BroadcastTable(table.high);
    destination0.low = _mm256_xor_si256(destination0.low,
        AVX2FF8ProductVector(source0.low, low_table, high_table));
    destination0.high = _mm256_xor_si256(destination0.high,
        AVX2FF8ProductVector(source0.high, low_table, high_table));
    destination1.low = _mm256_xor_si256(destination1.low,
        AVX2FF8ProductVector(source1.low, low_table, high_table));
    destination1.high = _mm256_xor_si256(destination1.high,
        AVX2FF8ProductVector(source1.high, low_table, high_table));
    destination2.low = _mm256_xor_si256(destination2.low,
        AVX2FF8ProductVector(source2.low, low_table, high_table));
    destination2.high = _mm256_xor_si256(destination2.high,
        AVX2FF8ProductVector(source2.high, low_table, high_table));
    destination3.low = _mm256_xor_si256(destination3.low,
        AVX2FF8ProductVector(source3.low, low_table, high_table));
    destination3.high = _mm256_xor_si256(destination3.high,
        AVX2FF8ProductVector(source3.high, low_table, high_table));
}

static LEO2_AVX2_FORCE_INLINE void AVX2FF8T8IFFTRadix4(
    AVX2FF8T8Pair& value0,
    AVX2FF8T8Pair& value1,
    AVX2FF8T8Pair& value2,
    AVX2FF8T8Pair& value3,
    uint16_t log01,
    uint16_t log23,
    uint16_t log02)
{
    static const uint16_t kZeroSkew = 255;
    AVX2FF8T8Xor(value1, value0);
    if (log01 != kZeroSkew)
        AVX2FF8T8MulAdd(value0, value1, log01);
    AVX2FF8T8Xor(value3, value2);
    if (log23 != kZeroSkew)
        AVX2FF8T8MulAdd(value2, value3, log23);
    AVX2FF8T8Xor(value2, value0);
    AVX2FF8T8Xor(value3, value1);
    if (log02 != kZeroSkew)
        AVX2FF8T8MulAdd2(value0, value2, value1, value3, log02);
}

static LEO2_AVX2_FORCE_INLINE void AVX2FF8T8IFFTDistance4(
    AVX2FF8T8Pair values[8],
    uint16_t log)
{
    static const uint16_t kZeroSkew = 255;
    AVX2FF8T8Xor(values[4], values[0]);
    AVX2FF8T8Xor(values[5], values[1]);
    AVX2FF8T8Xor(values[6], values[2]);
    AVX2FF8T8Xor(values[7], values[3]);
    if (log != kZeroSkew)
        AVX2FF8T8MulAdd4(
            values[0], values[4], values[1], values[5],
            values[2], values[6], values[3], values[7], log);
}

static LEO2_AVX2_FORCE_INLINE void AVX2FF8T8FFTRadix4Distance2(
    AVX2FF8T8Pair values[8],
    uint16_t log01,
    uint16_t log23,
    uint16_t log02)
{
    static const uint16_t kZeroSkew = 255;
    if (log02 != kZeroSkew)
        AVX2FF8T8MulAdd4(
            values[0], values[4], values[1], values[5],
            values[2], values[6], values[3], values[7], log02);
    AVX2FF8T8Xor(values[4], values[0]);
    AVX2FF8T8Xor(values[5], values[1]);
    AVX2FF8T8Xor(values[6], values[2]);
    AVX2FF8T8Xor(values[7], values[3]);
    if (log01 != kZeroSkew)
        AVX2FF8T8MulAdd2(
            values[0], values[2], values[1], values[3], log01);
    AVX2FF8T8Xor(values[2], values[0]);
    AVX2FF8T8Xor(values[3], values[1]);
    if (log23 != kZeroSkew)
        AVX2FF8T8MulAdd2(
            values[4], values[6], values[5], values[7], log23);
    AVX2FF8T8Xor(values[6], values[4]);
    AVX2FF8T8Xor(values[7], values[5]);
}

static LEO2_AVX2_FORCE_INLINE void AVX2FF8T8FFT2(
    AVX2FF8T8Pair& value0,
    AVX2FF8T8Pair& value1,
    uint16_t log)
{
    static const uint16_t kZeroSkew = 255;
    if (log != kZeroSkew)
        AVX2FF8T8MulAdd(value0, value1, log);
    AVX2FF8T8Xor(value1, value0);
}

static LEO2_AVX2_FORCE_INLINE void AVX2FF8T8ScalarIFFT2(
    uint8_t& value0,
    uint8_t& value1,
    uint16_t log)
{
    value1 ^= value0;
    if (log != 255)
        value0 ^= FF8Product(log, value1);
}

static LEO2_AVX2_FORCE_INLINE void AVX2FF8T8ScalarFFT2(
    uint8_t& value0,
    uint8_t& value1,
    uint16_t log)
{
    if (log != 255)
        value0 ^= FF8Product(log, value1);
    value1 ^= value0;
}

#if defined(_MSC_VER)
# define LEO2_AVX2_T8_NOINLINE __declspec(noinline)
#elif defined(__GNUC__) || defined(__clang__)
# define LEO2_AVX2_T8_NOINLINE __attribute__((noinline))
#else
# define LEO2_AVX2_T8_NOINLINE
#endif

static LEO2_AVX2_T8_NOINLINE void AVX2FF8HighEncodeT8ScalarTail(
    const uint8_t* const input[8],
    uint8_t* const output[8],
    uint64_t lane7_offset_mask,
    const uint8_t* inverse_skew,
    const uint8_t* forward_skew,
    uint64_t offset,
    uint64_t byte_count)
{
    while (offset < byte_count)
    {
        uint8_t values[8];
        for (unsigned lane = 0; lane < 7; ++lane)
            values[lane] = input[lane][offset];
        values[7] = input[7][offset & lane7_offset_mask];

        AVX2FF8T8ScalarIFFT2(values[0], values[1], inverse_skew[1]);
        AVX2FF8T8ScalarIFFT2(values[2], values[3], inverse_skew[3]);
        AVX2FF8T8ScalarIFFT2(values[0], values[2], inverse_skew[2]);
        AVX2FF8T8ScalarIFFT2(values[1], values[3], inverse_skew[2]);
        AVX2FF8T8ScalarIFFT2(values[4], values[5], inverse_skew[5]);
        AVX2FF8T8ScalarIFFT2(values[6], values[7], inverse_skew[7]);
        AVX2FF8T8ScalarIFFT2(values[4], values[6], inverse_skew[6]);
        AVX2FF8T8ScalarIFFT2(values[5], values[7], inverse_skew[6]);
        for (unsigned lane = 0; lane < 4; ++lane)
            AVX2FF8T8ScalarIFFT2(
                values[lane], values[lane + 4], inverse_skew[4]);

        AVX2FF8T8ScalarFFT2(values[0], values[4], forward_skew[4]);
        AVX2FF8T8ScalarFFT2(values[1], values[5], forward_skew[4]);
        AVX2FF8T8ScalarFFT2(values[2], values[6], forward_skew[4]);
        AVX2FF8T8ScalarFFT2(values[3], values[7], forward_skew[4]);
        AVX2FF8T8ScalarFFT2(values[0], values[2], forward_skew[2]);
        AVX2FF8T8ScalarFFT2(values[1], values[3], forward_skew[2]);
        AVX2FF8T8ScalarFFT2(values[4], values[6], forward_skew[6]);
        AVX2FF8T8ScalarFFT2(values[5], values[7], forward_skew[6]);
        AVX2FF8T8ScalarFFT2(values[0], values[1], forward_skew[1]);
        AVX2FF8T8ScalarFFT2(values[2], values[3], forward_skew[3]);
        AVX2FF8T8ScalarFFT2(values[4], values[5], forward_skew[5]);
        AVX2FF8T8ScalarFFT2(values[6], values[7], forward_skew[7]);

        for (unsigned lane = 0; lane < 8; ++lane)
            output[lane][offset] = values[lane];
        ++offset;
    }
}

#undef LEO2_AVX2_T8_NOINLINE

static void AVX2FF8HighEncodeT8(
    const void* const* data,
    void* const* work,
    bool shortened,
    const uint8_t* inverse_skew,
    const uint8_t* forward_skew,
    uint64_t byte_count)
{
    alignas(32) static const uint8_t kZeroShard[64] = {};
    const uint8_t* input[8];
    uint8_t* output[8];
    for (unsigned lane = 0; lane < 7; ++lane)
    {
        input[lane] = static_cast<const uint8_t*>(data[lane]);
        output[lane] = static_cast<uint8_t*>(work[lane]);
    }
    input[7] = shortened
        ? kZeroShard : static_cast<const uint8_t*>(data[7]);
    output[7] = static_cast<uint8_t*>(work[7]);
    const uint64_t lane7_offset_mask = shortened ? 0 : ~uint64_t(0);

    uint64_t offset = 0;
    while (byte_count - offset >= 64)
    {
        AVX2FF8T8Pair values[8];
        for (unsigned lane = 0; lane < 7; ++lane)
        {
            values[lane].low = _mm256_loadu_si256(
                reinterpret_cast<const __m256i*>(input[lane] + offset));
            values[lane].high = _mm256_loadu_si256(
                reinterpret_cast<const __m256i*>(
                    input[lane] + offset + 32));
        }
        const uint64_t lane7_offset = offset & lane7_offset_mask;
        values[7].low = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(input[7] + lane7_offset));
        values[7].high = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(
                input[7] + lane7_offset + 32));

        AVX2FF8T8IFFTRadix4(
            values[0], values[1], values[2], values[3],
            inverse_skew[1], inverse_skew[3], inverse_skew[2]);
        AVX2FF8T8IFFTRadix4(
            values[4], values[5], values[6], values[7],
            inverse_skew[5], inverse_skew[7], inverse_skew[6]);
        AVX2FF8T8IFFTDistance4(values, inverse_skew[4]);

        AVX2FF8T8FFTRadix4Distance2(values,
            forward_skew[2], forward_skew[6], forward_skew[4]);
        AVX2FF8T8FFT2(values[0], values[1], forward_skew[1]);
        AVX2FF8T8FFT2(values[2], values[3], forward_skew[3]);
        AVX2FF8T8FFT2(values[4], values[5], forward_skew[5]);
        AVX2FF8T8FFT2(values[6], values[7], forward_skew[7]);

        for (unsigned lane = 0; lane < 8; ++lane)
        {
            _mm256_storeu_si256(
                reinterpret_cast<__m256i*>(output[lane] + offset),
                values[lane].low);
            _mm256_storeu_si256(
                reinterpret_cast<__m256i*>(output[lane] + offset + 32),
                values[lane].high);
        }
        offset += 64;
    }

    if (offset < byte_count)
        AVX2FF8HighEncodeT8ScalarTail(
            input, output, lane7_offset_mask,
            inverse_skew, forward_skew, offset, byte_count);
}

#endif // LEO2_AVX512_VARIANT

template<bool Inverse, bool AllNonzero>
static void AVX2FF8Butterfly4RangePreparedImpl(
    void* const* work,
    uint32_t distance,
    uint16_t log01,
    uint16_t log23,
    uint16_t log02,
    uint64_t byte_count)
{
    static const uint16_t kZeroSkew = 255;
    __m256i low01 = _mm256_setzero_si256();
    __m256i high01 = _mm256_setzero_si256();
    __m256i low23 = _mm256_setzero_si256();
    __m256i high23 = _mm256_setzero_si256();
    __m256i low02 = _mm256_setzero_si256();
    __m256i high02 = _mm256_setzero_si256();
    if (AllNonzero || log01 != kZeroSkew)
    {
        low01 = BroadcastTable(FF8Tables[log01].low);
        high01 = BroadcastTable(FF8Tables[log01].high);
    }
    if (AllNonzero || log23 != kZeroSkew)
    {
        low23 = BroadcastTable(FF8Tables[log23].low);
        high23 = BroadcastTable(FF8Tables[log23].high);
    }
    if (AllNonzero || log02 != kZeroSkew)
    {
        low02 = BroadcastTable(FF8Tables[log02].low);
        high02 = BroadcastTable(FF8Tables[log02].high);
    }

    for (uint32_t lane = 0; lane < distance; ++lane)
    {
        uint8_t* value0 = static_cast<uint8_t*>(work[lane]);
        uint8_t* value1 = static_cast<uint8_t*>(work[lane + distance]);
        uint8_t* value2 = static_cast<uint8_t*>(work[lane + distance * 2U]);
        uint8_t* value3 = static_cast<uint8_t*>(work[lane + distance * 3U]);
        uint64_t remaining = byte_count;
        while (remaining >= 32)
        {
            __m256i x0 = _mm256_loadu_si256(
                reinterpret_cast<const __m256i*>(value0));
            __m256i x1 = _mm256_loadu_si256(
                reinterpret_cast<const __m256i*>(value1));
            __m256i x2 = _mm256_loadu_si256(
                reinterpret_cast<const __m256i*>(value2));
            __m256i x3 = _mm256_loadu_si256(
                reinterpret_cast<const __m256i*>(value3));
            if (Inverse)
            {
                x1 = _mm256_xor_si256(x1, x0);
                if (AllNonzero || log01 != kZeroSkew)
                    x0 = _mm256_xor_si256(x0,
                        AVX2FF8ProductVector(x1, low01, high01));
                x3 = _mm256_xor_si256(x3, x2);
                if (AllNonzero || log23 != kZeroSkew)
                    x2 = _mm256_xor_si256(x2,
                        AVX2FF8ProductVector(x3, low23, high23));
                x2 = _mm256_xor_si256(x2, x0);
                x3 = _mm256_xor_si256(x3, x1);
                if (AllNonzero || log02 != kZeroSkew)
                {
                    x0 = _mm256_xor_si256(x0,
                        AVX2FF8ProductVector(x2, low02, high02));
                    x1 = _mm256_xor_si256(x1,
                        AVX2FF8ProductVector(x3, low02, high02));
                }
            }
            else
            {
                if (AllNonzero || log02 != kZeroSkew)
                {
                    x0 = _mm256_xor_si256(x0,
                        AVX2FF8ProductVector(x2, low02, high02));
                    x1 = _mm256_xor_si256(x1,
                        AVX2FF8ProductVector(x3, low02, high02));
                }
                x2 = _mm256_xor_si256(x2, x0);
                x3 = _mm256_xor_si256(x3, x1);
                if (AllNonzero || log01 != kZeroSkew)
                    x0 = _mm256_xor_si256(x0,
                        AVX2FF8ProductVector(x1, low01, high01));
                x1 = _mm256_xor_si256(x1, x0);
                if (AllNonzero || log23 != kZeroSkew)
                    x2 = _mm256_xor_si256(x2,
                        AVX2FF8ProductVector(x3, low23, high23));
                x3 = _mm256_xor_si256(x3, x2);
            }
            _mm256_storeu_si256(reinterpret_cast<__m256i*>(value0), x0);
            _mm256_storeu_si256(reinterpret_cast<__m256i*>(value1), x1);
            _mm256_storeu_si256(reinterpret_cast<__m256i*>(value2), x2);
            _mm256_storeu_si256(reinterpret_cast<__m256i*>(value3), x3);
            value0 += 32;
            value1 += 32;
            value2 += 32;
            value3 += 32;
            remaining -= 32;
        }
        while (remaining-- != 0)
        {
            uint8_t x0 = *value0;
            uint8_t x1 = *value1;
            uint8_t x2 = *value2;
            uint8_t x3 = *value3;
            if (Inverse)
            {
                x1 ^= x0;
                if (AllNonzero || log01 != kZeroSkew)
                    x0 ^= FF8Product(log01, x1);
                x3 ^= x2;
                if (AllNonzero || log23 != kZeroSkew)
                    x2 ^= FF8Product(log23, x3);
                x2 ^= x0;
                x3 ^= x1;
                if (AllNonzero || log02 != kZeroSkew)
                {
                    x0 ^= FF8Product(log02, x2);
                    x1 ^= FF8Product(log02, x3);
                }
            }
            else
            {
                if (AllNonzero || log02 != kZeroSkew)
                {
                    x0 ^= FF8Product(log02, x2);
                    x1 ^= FF8Product(log02, x3);
                }
                x2 ^= x0;
                x3 ^= x1;
                if (AllNonzero || log01 != kZeroSkew)
                    x0 ^= FF8Product(log01, x1);
                x1 ^= x0;
                if (AllNonzero || log23 != kZeroSkew)
                    x2 ^= FF8Product(log23, x3);
                x3 ^= x2;
            }
            *value0++ = x0;
            *value1++ = x1;
            *value2++ = x2;
            *value3++ = x3;
        }
    }
}

template<bool Inverse>
static void AVX2FF8Butterfly4RangePrepared(
    void* const* work,
    uint32_t distance,
    uint16_t log01,
    uint16_t log23,
    uint16_t log02,
    uint64_t byte_count)
{
    static const uint16_t kZeroSkew = 255;
    if (log01 != kZeroSkew && log23 != kZeroSkew && log02 != kZeroSkew)
        AVX2FF8Butterfly4RangePreparedImpl<Inverse, true>(
            work, distance, log01, log23, log02, byte_count);
    else
        AVX2FF8Butterfly4RangePreparedImpl<Inverse, false>(
            work, distance, log01, log23, log02, byte_count);
}

#if defined(LEO_HAS_FF8) && !defined(LEO2_AVX512_VARIANT) && \
    !defined(LEO2_GFNI_VARIANT)

#if defined(_MSC_VER)
#define LEO2_T32_FINAL_NOINLINE __declspec(noinline)
#elif (defined(__GNUC__) || defined(__clang__)) && !defined(__APPLE__)
#define LEO2_T32_FINAL_NOINLINE \
    __attribute__((noinline, section(".text.leo2_t32_final")))
#elif defined(__GNUC__) || defined(__clang__)
#define LEO2_T32_FINAL_NOINLINE __attribute__((noinline))
#else
#define LEO2_T32_FINAL_NOINLINE
#endif

// Exact 64-byte T=32 final inverse layer.  The low/high tables remain live
// across every pair instead of being reloaded by sixteen leaf callbacks.
static LEO2_T32_FINAL_NOINLINE void AVX2FF8IFFTButterfly2Range(
    void* const* work,
    void* const* xor_output,
    unsigned distance,
    uint16_t multiplier_log,
    uint64_t byte_count)
{
    static const uint16_t kZeroSkew = 255;
    LEO_DEBUG_ASSERT(work != NULL);
    LEO_DEBUG_ASSERT(distance != 0);
    LEO_DEBUG_ASSERT(multiplier_log != kZeroSkew);
    LEO_DEBUG_ASSERT(byte_count == 64);
    if (!work || distance == 0 || multiplier_log == kZeroSkew ||
        byte_count != 64)
        return;

    const FF8NibbleTable& table = FF8Tables[multiplier_log];
    const __m256i low_table = BroadcastTable(table.low);
    const __m256i high_table = BroadcastTable(table.high);

    if (xor_output)
    {
        for (unsigned lane = 0; lane < distance; ++lane)
        {
            const uint8_t* x = static_cast<const uint8_t*>(work[lane]);
            const uint8_t* y = static_cast<const uint8_t*>(
                work[lane + distance]);
            uint8_t* out_x = static_cast<uint8_t*>(xor_output[lane]);
            uint8_t* out_y = static_cast<uint8_t*>(
                xor_output[lane + distance]);

            const __m256i x0 = _mm256_loadu_si256(
                reinterpret_cast<const __m256i*>(x));
            const __m256i x1 = _mm256_loadu_si256(
                reinterpret_cast<const __m256i*>(x + 32));
            const __m256i y0 = _mm256_xor_si256(x0,
                _mm256_loadu_si256(reinterpret_cast<const __m256i*>(y)));
            const __m256i y1 = _mm256_xor_si256(x1,
                _mm256_loadu_si256(
                    reinterpret_cast<const __m256i*>(y + 32)));
            const __m256i transformed_x0 = _mm256_xor_si256(x0,
                AVX2FF8ProductVector(y0, low_table, high_table));
            const __m256i transformed_x1 = _mm256_xor_si256(x1,
                AVX2FF8ProductVector(y1, low_table, high_table));
            const __m256i result_x0 = _mm256_xor_si256(transformed_x0,
                _mm256_loadu_si256(reinterpret_cast<const __m256i*>(out_x)));
            const __m256i result_x1 = _mm256_xor_si256(transformed_x1,
                _mm256_loadu_si256(
                    reinterpret_cast<const __m256i*>(out_x + 32)));
            const __m256i result_y0 = _mm256_xor_si256(y0,
                _mm256_loadu_si256(reinterpret_cast<const __m256i*>(out_y)));
            const __m256i result_y1 = _mm256_xor_si256(y1,
                _mm256_loadu_si256(
                    reinterpret_cast<const __m256i*>(out_y + 32)));
            _mm256_storeu_si256(
                reinterpret_cast<__m256i*>(out_x), result_x0);
            _mm256_storeu_si256(
                reinterpret_cast<__m256i*>(out_x + 32), result_x1);
            _mm256_storeu_si256(
                reinterpret_cast<__m256i*>(out_y), result_y0);
            _mm256_storeu_si256(
                reinterpret_cast<__m256i*>(out_y + 32), result_y1);
        }
    }
    else
    {
        for (unsigned lane = 0; lane < distance; ++lane)
        {
            uint8_t* x = static_cast<uint8_t*>(work[lane]);
            uint8_t* y = static_cast<uint8_t*>(work[lane + distance]);
            __m256i x0 = _mm256_loadu_si256(
                reinterpret_cast<const __m256i*>(x));
            __m256i x1 = _mm256_loadu_si256(
                reinterpret_cast<const __m256i*>(x + 32));
            __m256i y0 = _mm256_xor_si256(x0,
                _mm256_loadu_si256(reinterpret_cast<const __m256i*>(y)));
            __m256i y1 = _mm256_xor_si256(x1,
                _mm256_loadu_si256(
                    reinterpret_cast<const __m256i*>(y + 32)));
            x0 = _mm256_xor_si256(x0,
                AVX2FF8ProductVector(y0, low_table, high_table));
            x1 = _mm256_xor_si256(x1,
                AVX2FF8ProductVector(y1, low_table, high_table));
            _mm256_storeu_si256(reinterpret_cast<__m256i*>(x), x0);
            _mm256_storeu_si256(
                reinterpret_cast<__m256i*>(x + 32), x1);
            _mm256_storeu_si256(reinterpret_cast<__m256i*>(y), y0);
            _mm256_storeu_si256(
                reinterpret_cast<__m256i*>(y + 32), y1);
        }
    }
}

#undef LEO2_T32_FINAL_NOINLINE
#endif

#if defined(LEO2_AVX512_VARIANT)

static void AVX2FF8IFFTButterfly2RangePrepared(
    void* const* work,
    uint32_t distance,
    uint16_t log,
    uint64_t byte_count)
{
    static const uint16_t kZeroSkew = 255;
    if (log == kZeroSkew)
    {
        for (uint32_t lane = 0; lane < distance; ++lane)
            AVX2XorMemory(
                work[lane + distance], work[lane], byte_count);
        return;
    }
    const __m256i low = BroadcastTable(FF8Tables[log].low);
    const __m256i high = BroadcastTable(FF8Tables[log].high);
    for (uint32_t lane = 0; lane < distance; ++lane)
    {
        uint8_t* x = static_cast<uint8_t*>(work[lane]);
        uint8_t* y = static_cast<uint8_t*>(work[lane + distance]);
        uint64_t remaining = byte_count;
        while (remaining >= 32)
        {
            __m256i x_value = _mm256_loadu_si256(
                reinterpret_cast<const __m256i*>(x));
            __m256i y_value = _mm256_loadu_si256(
                reinterpret_cast<const __m256i*>(y));
            y_value = _mm256_xor_si256(y_value, x_value);
            x_value = _mm256_xor_si256(x_value,
                AVX2FF8ProductVector(y_value, low, high));
            _mm256_storeu_si256(reinterpret_cast<__m256i*>(x), x_value);
            _mm256_storeu_si256(reinterpret_cast<__m256i*>(y), y_value);
            x += 32;
            y += 32;
            remaining -= 32;
        }
        while (remaining-- != 0)
        {
            *y ^= *x;
            *x ^= FF8Product(log, *y);
            ++x;
            ++y;
        }
    }
}

static void AVX2FF8HighEncodeOneBlock(
    const void* const* data,
    void* const* work,
    uint32_t side_and_flags,
    const uint8_t* inverse_skew,
    const uint8_t* forward_skew,
    uint64_t byte_count)
{
    static const uint16_t kZeroSkew = 255;
    const bool shortened =
        (side_and_flags & kFF8HighEncodeShortenedInput) != 0;
    const uint32_t side =
        side_and_flags & ~kFF8HighEncodeShortenedInput;
    const uint32_t data_count = side - static_cast<uint32_t>(shortened);
    LEO_DEBUG_ASSERT(side == 8 || side == 16 || side == 32 || side == 64);

    if (side == 8)
    {
        AVX2FF8HighEncodeT8(
            data, work, shortened, inverse_skew, forward_skew, byte_count);
        return;
    }

    // Consume the systematic shards directly through the first two inverse
    // layers, exactly as IFFT_DIT_Encoder does for a complete input block.  A
    // one-coordinate shortened tail changes only the final four-row boundary;
    // the weighted primitive supplies mathematical zero without reading a
    // nonexistent public pointer and leaves the remaining traversal regular.
    for (uint32_t r = 0; r < side; r += 4)
    {
        if (!shortened || r + 4U <= data_count)
        {
            AVX2FF8IFFTButterfly4Out(
                data[r], data[r + 1], data[r + 2], data[r + 3],
                work[r], work[r + 1], work[r + 2], work[r + 3],
                inverse_skew[r + 1], inverse_skew[r + 3],
                inverse_skew[r + 2], byte_count);
        }
        else
        {
            LEO_DEBUG_ASSERT(r + 3U == data_count);
            AVX2FF8IFFTButterfly4OutLastZero(
                data[r], data[r + 1], data[r + 2],
                work[r], work[r + 1], work[r + 2], work[r + 3],
                inverse_skew[r + 1], inverse_skew[r + 3],
                inverse_skew[r + 2], byte_count);
        }
    }

    uint32_t distance = 4;
    uint32_t distance4 = 16;
    for (; distance4 <= side;
         distance = distance4, distance4 <<= 2)
    {
        for (uint32_t r = 0; r < side; r += distance4)
        {
            const uint32_t end = r + distance;
            AVX2FF8Butterfly4RangePrepared<true>(
                work + r, distance,
                inverse_skew[end], inverse_skew[end + distance * 2U],
                inverse_skew[end + distance], byte_count);
        }
    }
    if (distance < side)
    {
        AVX2FF8IFFTButterfly2RangePrepared(
            work, distance, inverse_skew[distance], byte_count);
    }

    // Evaluate the coefficient block on the legacy parity coset.  The
    // promoted production schedule uses the same fused radix-four operation
    // at every complete stage for these sizes.
    distance4 = side;
    distance = side >> 2;
    for (; distance != 0;
         distance4 = distance, distance >>= 2)
    {
        for (uint32_t r = 0; r < side; r += distance4)
        {
            const uint32_t end = r + distance;
            AVX2FF8Butterfly4RangePrepared<false>(
                work + r, distance,
                forward_skew[end], forward_skew[end + distance * 2U],
                forward_skew[end + distance], byte_count);
        }
    }
    if (distance4 == 2)
    {
        for (uint32_t r = 0; r < side; r += 2)
        {
            const uint16_t log = forward_skew[r + 1];
            if (log == kZeroSkew)
                AVX2XorMemory(work[r + 1], work[r], byte_count);
            else
                AVX2FF8Butterfly2<false>(
                    work[r], work[r + 1], log, byte_count);
        }
    }
}

#endif // LEO2_AVX512_VARIANT

#endif // LEO_HAS_FF8

#ifdef LEO_HAS_FF16
static LEO2_AVX2_FORCE_INLINE void AVX2FF16MultiplyAddPair(
    __m256i& destination_low,
    __m256i& destination_high,
    __m256i source_low,
    __m256i source_high,
    const FF16Table& table)
{
#if defined(LEO2_GFNI_VARIANT)
    const __m256i low_tables[4] = {
        BroadcastAffine(table.block[0]), BroadcastAffine(table.block[1]),
        BroadcastAffine(table.block[2]), BroadcastAffine(table.block[3])
    };
    // ProductVectors ignores its high_tables operand in this variant.
    const __m256i* const high_tables = low_tables;
#else
    const __m256i low_tables[4] = {
        BroadcastTable(table.low[0]), BroadcastTable(table.low[1]),
        BroadcastTable(table.low[2]), BroadcastTable(table.low[3])
    };
    const __m256i high_tables[4] = {
        BroadcastTable(table.high[0]), BroadcastTable(table.high[1]),
        BroadcastTable(table.high[2]), BroadcastTable(table.high[3])
    };
#endif
    __m256i product_low;
    __m256i product_high;
    AVX2FF16ProductVectors(source_low, source_high,
        low_tables, high_tables, product_low, product_high);
    destination_low = _mm256_xor_si256(destination_low, product_low);
    destination_high = _mm256_xor_si256(destination_high, product_high);
}

template<bool Inverse>
static void AVX2FF16Butterfly4Split(
    void* value0, void* value1, void* value2, void* value3,
    uint16_t log01, uint16_t log23, uint16_t log02,
    uint64_t byte_count)
{
    static const uint16_t kZeroSkew = 65535;
    if (Inverse)
    {
        if (log01 == kZeroSkew)
            AVX2XorMemory(value1, value0, byte_count);
        else
            AVX2FF16Butterfly2<true>(value0, value1, log01, byte_count);
        if (log23 == kZeroSkew)
            AVX2XorMemory(value3, value2, byte_count);
        else
            AVX2FF16Butterfly2<true>(value2, value3, log23, byte_count);
        if (log02 == kZeroSkew)
        {
            AVX2XorMemory(value2, value0, byte_count);
            AVX2XorMemory(value3, value1, byte_count);
        }
        else
        {
            AVX2FF16Butterfly2<true>(value0, value2, log02, byte_count);
            AVX2FF16Butterfly2<true>(value1, value3, log02, byte_count);
        }
    }
    else
    {
        if (log02 == kZeroSkew)
        {
            AVX2XorMemory(value2, value0, byte_count);
            AVX2XorMemory(value3, value1, byte_count);
        }
        else
        {
            AVX2FF16Butterfly2<false>(value0, value2, log02, byte_count);
            AVX2FF16Butterfly2<false>(value1, value3, log02, byte_count);
        }
        if (log01 == kZeroSkew)
            AVX2XorMemory(value1, value0, byte_count);
        else
            AVX2FF16Butterfly2<false>(value0, value1, log01, byte_count);
        if (log23 == kZeroSkew)
            AVX2XorMemory(value3, value2, byte_count);
        else
            AVX2FF16Butterfly2<false>(value2, value3, log23, byte_count);
    }
}

template<bool Inverse>
static void AVX2FF16Butterfly4(
    void* value0, void* value1, void* value2, void* value3,
    uint16_t log01, uint16_t log23, uint16_t log02,
    uint64_t byte_count)
{
    static const uint16_t kZeroSkew = 65535;
    // Keep small working sets fused.  This also covers public tails such as
    // 66 bytes, which staging rounds to 128 bytes.  Larger shards use the
    // split schedule to bound register pressure and table residency.
    //
    // The GFNI variant has no such bound.  A nibble multiplier needs eight
    // broadcast table vectors, so three of them plus eight data vectors do not
    // fit in sixteen ymm; an affine multiplier needs four, they are
    // re-broadcast at each use rather than held, and the live set is eight
    // data vectors plus four matrices plus two products.  Fusing at every size
    // loads and stores the four shards once per radix-four group instead of
    // four times, which measured 1.72x to 1.89x against the split schedule
    // from 1 KiB to 256 KiB per shard.
#if defined(LEO2_GFNI_VARIANT)
    if (byte_count < 64)
#else
    if (byte_count < 64 || byte_count > 128)
#endif
    {
        AVX2FF16Butterfly4Split<Inverse>(
            value0, value1, value2, value3,
            log01, log23, log02, byte_count);
        return;
    }
    uint8_t* values[4] = {
        static_cast<uint8_t*>(value0), static_cast<uint8_t*>(value1),
        static_cast<uint8_t*>(value2), static_cast<uint8_t*>(value3)
    };
    uint64_t offset = 0;
    while (byte_count - offset >= 64)
    {
        __m256i low[4];
        __m256i high[4];
        for (unsigned lane = 0; lane < 4; ++lane)
        {
            low[lane] = _mm256_loadu_si256(
                reinterpret_cast<const __m256i*>(values[lane] + offset));
            high[lane] = _mm256_loadu_si256(
                reinterpret_cast<const __m256i*>(
                    values[lane] + offset + 32));
        }

        if (Inverse)
        {
            low[1] = _mm256_xor_si256(low[1], low[0]);
            high[1] = _mm256_xor_si256(high[1], high[0]);
            if (log01 != kZeroSkew)
                AVX2FF16MultiplyAddPair(low[0], high[0], low[1], high[1],
                    FF16Tables[log01]);

            low[3] = _mm256_xor_si256(low[3], low[2]);
            high[3] = _mm256_xor_si256(high[3], high[2]);
            if (log23 != kZeroSkew)
                AVX2FF16MultiplyAddPair(low[2], high[2], low[3], high[3],
                    FF16Tables[log23]);

            low[2] = _mm256_xor_si256(low[2], low[0]);
            high[2] = _mm256_xor_si256(high[2], high[0]);
            low[3] = _mm256_xor_si256(low[3], low[1]);
            high[3] = _mm256_xor_si256(high[3], high[1]);
            if (log02 != kZeroSkew)
            {
                AVX2FF16MultiplyAddPair(low[0], high[0], low[2], high[2],
                    FF16Tables[log02]);
                AVX2FF16MultiplyAddPair(low[1], high[1], low[3], high[3],
                    FF16Tables[log02]);
            }
        }
        else
        {
            if (log02 != kZeroSkew)
            {
                AVX2FF16MultiplyAddPair(low[0], high[0], low[2], high[2],
                    FF16Tables[log02]);
                AVX2FF16MultiplyAddPair(low[1], high[1], low[3], high[3],
                    FF16Tables[log02]);
            }
            low[2] = _mm256_xor_si256(low[2], low[0]);
            high[2] = _mm256_xor_si256(high[2], high[0]);
            low[3] = _mm256_xor_si256(low[3], low[1]);
            high[3] = _mm256_xor_si256(high[3], high[1]);

            if (log01 != kZeroSkew)
                AVX2FF16MultiplyAddPair(low[0], high[0], low[1], high[1],
                    FF16Tables[log01]);
            low[1] = _mm256_xor_si256(low[1], low[0]);
            high[1] = _mm256_xor_si256(high[1], high[0]);

            if (log23 != kZeroSkew)
                AVX2FF16MultiplyAddPair(low[2], high[2], low[3], high[3],
                    FF16Tables[log23]);
            low[3] = _mm256_xor_si256(low[3], low[2]);
            high[3] = _mm256_xor_si256(high[3], high[2]);
        }

        for (unsigned lane = 0; lane < 4; ++lane)
        {
            _mm256_storeu_si256(
                reinterpret_cast<__m256i*>(values[lane] + offset), low[lane]);
            _mm256_storeu_si256(reinterpret_cast<__m256i*>(
                values[lane] + offset + 32), high[lane]);
        }
        offset += 64;
    }

    const uint64_t residual = byte_count - offset;
    if (residual == 0)
        return;
    AVX2FF16Butterfly4Split<Inverse>(
        values[0] + offset, values[1] + offset,
        values[2] + offset, values[3] + offset,
        log01, log23, log02, residual);
}

static void AVX2FF16IFFTButterfly4(
    void* value0, void* value1, void* value2, void* value3,
    uint16_t log01, uint16_t log23, uint16_t log02,
    uint64_t byte_count)
{
    AVX2FF16Butterfly4<true>(value0, value1, value2, value3,
        log01, log23, log02, byte_count);
}

#undef LEO2_AVX2_FORCE_INLINE

static void AVX2FF16FFTButterfly4(
    void* value0, void* value1, void* value2, void* value3,
    uint16_t log01, uint16_t log23, uint16_t log02,
    uint64_t byte_count)
{
    AVX2FF16Butterfly4<false>(value0, value1, value2, value3,
        log01, log23, log02, byte_count);
}

template<bool Inverse>
static void AVX2FF16Butterfly4Out(
    const void* input0_pointer, const void* input1_pointer,
    const void* input2_pointer, const void* input3_pointer,
    void* output0_pointer, void* output1_pointer,
    void* output2_pointer, void* output3_pointer,
    uint16_t log01, uint16_t log23, uint16_t log02,
    uint64_t byte_count)
{
    static const uint16_t kZeroSkew = 65535;
    const uint8_t* inputs[4] = {
        static_cast<const uint8_t*>(input0_pointer),
        static_cast<const uint8_t*>(input1_pointer),
        static_cast<const uint8_t*>(input2_pointer),
        static_cast<const uint8_t*>(input3_pointer)
    };
    uint8_t* outputs[4] = {
        static_cast<uint8_t*>(output0_pointer),
        static_cast<uint8_t*>(output1_pointer),
        static_cast<uint8_t*>(output2_pointer),
        static_cast<uint8_t*>(output3_pointer)
    };
    uint64_t offset = 0;
    while (byte_count - offset >= 64)
    {
        __m256i low[4];
        __m256i high[4];
        for (unsigned lane = 0; lane < 4; ++lane)
        {
            low[lane] = _mm256_loadu_si256(
                reinterpret_cast<const __m256i*>(inputs[lane] + offset));
            high[lane] = _mm256_loadu_si256(
                reinterpret_cast<const __m256i*>(
                    inputs[lane] + offset + 32));
        }
        if (Inverse)
        {
            low[1] = _mm256_xor_si256(low[1], low[0]);
            high[1] = _mm256_xor_si256(high[1], high[0]);
            if (log01 != kZeroSkew)
                AVX2FF16MultiplyAddPair(low[0], high[0],
                    low[1], high[1], FF16Tables[log01]);
            low[3] = _mm256_xor_si256(low[3], low[2]);
            high[3] = _mm256_xor_si256(high[3], high[2]);
            if (log23 != kZeroSkew)
                AVX2FF16MultiplyAddPair(low[2], high[2],
                    low[3], high[3], FF16Tables[log23]);
            low[2] = _mm256_xor_si256(low[2], low[0]);
            high[2] = _mm256_xor_si256(high[2], high[0]);
            low[3] = _mm256_xor_si256(low[3], low[1]);
            high[3] = _mm256_xor_si256(high[3], high[1]);
            if (log02 != kZeroSkew)
            {
                AVX2FF16MultiplyAddPair(low[0], high[0],
                    low[2], high[2], FF16Tables[log02]);
                AVX2FF16MultiplyAddPair(low[1], high[1],
                    low[3], high[3], FF16Tables[log02]);
            }
        }
        else
        {
            if (log02 != kZeroSkew)
            {
                AVX2FF16MultiplyAddPair(low[0], high[0],
                    low[2], high[2], FF16Tables[log02]);
                AVX2FF16MultiplyAddPair(low[1], high[1],
                    low[3], high[3], FF16Tables[log02]);
            }
            low[2] = _mm256_xor_si256(low[2], low[0]);
            high[2] = _mm256_xor_si256(high[2], high[0]);
            low[3] = _mm256_xor_si256(low[3], low[1]);
            high[3] = _mm256_xor_si256(high[3], high[1]);
            if (log01 != kZeroSkew)
                AVX2FF16MultiplyAddPair(low[0], high[0],
                    low[1], high[1], FF16Tables[log01]);
            low[1] = _mm256_xor_si256(low[1], low[0]);
            high[1] = _mm256_xor_si256(high[1], high[0]);
            if (log23 != kZeroSkew)
                AVX2FF16MultiplyAddPair(low[2], high[2],
                    low[3], high[3], FF16Tables[log23]);
            low[3] = _mm256_xor_si256(low[3], low[2]);
            high[3] = _mm256_xor_si256(high[3], high[2]);
        }
        for (unsigned lane = 0; lane < 4; ++lane)
        {
            _mm256_storeu_si256(
                reinterpret_cast<__m256i*>(outputs[lane] + offset), low[lane]);
            _mm256_storeu_si256(reinterpret_cast<__m256i*>(
                outputs[lane] + offset + 32), high[lane]);
        }
        offset += 64;
    }
    const uint64_t symbols = (byte_count - offset) / 2;
    for (uint64_t i = 0; i < symbols; ++i)
    {
        uint16_t x[4];
        for (unsigned lane = 0; lane < 4; ++lane)
            x[lane] = static_cast<uint16_t>(inputs[lane][offset + i] |
                (static_cast<unsigned>(inputs[lane][offset + symbols + i]) << 8));
        if (Inverse)
        {
            x[1] ^= x[0];
            if (log01 != kZeroSkew)
                x[0] ^= FF16Product(log01, x[1]);
            x[3] ^= x[2];
            if (log23 != kZeroSkew)
                x[2] ^= FF16Product(log23, x[3]);
            x[2] ^= x[0];
            x[3] ^= x[1];
            if (log02 != kZeroSkew)
            {
                x[0] ^= FF16Product(log02, x[2]);
                x[1] ^= FF16Product(log02, x[3]);
            }
        }
        else
        {
            if (log02 != kZeroSkew)
            {
                x[0] ^= FF16Product(log02, x[2]);
                x[1] ^= FF16Product(log02, x[3]);
            }
            x[2] ^= x[0];
            x[3] ^= x[1];
            if (log01 != kZeroSkew)
                x[0] ^= FF16Product(log01, x[1]);
            x[1] ^= x[0];
            if (log23 != kZeroSkew)
                x[2] ^= FF16Product(log23, x[3]);
            x[3] ^= x[2];
        }
        for (unsigned lane = 0; lane < 4; ++lane)
        {
            outputs[lane][offset + i] = static_cast<uint8_t>(x[lane]);
            outputs[lane][offset + symbols + i] =
                static_cast<uint8_t>(x[lane] >> 8);
        }
    }
}

static void AVX2FF16IFFTButterfly4Out(
    const void* input0, const void* input1,
    const void* input2, const void* input3,
    void* output0, void* output1, void* output2, void* output3,
    uint16_t log01, uint16_t log23, uint16_t log02,
    uint64_t byte_count)
{
    AVX2FF16Butterfly4Out<true>(
        input0, input1, input2, input3,
        output0, output1, output2, output3,
        log01, log23, log02, byte_count);
}

static void AVX2FF16FFTButterfly4Out(
    const void* input0, const void* input1,
    const void* input2, const void* input3,
    void* output0, void* output1, void* output2, void* output3,
    uint16_t log01, uint16_t log23, uint16_t log02,
    uint64_t byte_count)
{
    AVX2FF16Butterfly4Out<false>(
        input0, input1, input2, input3,
        output0, output1, output2, output3,
        log01, log23, log02, byte_count);
}

#if !defined(LEO2_AVX512_VARIANT)
template<bool Inverse>
static void AVX2FF16Butterfly2RangePrepared(
    void* const* work,
    unsigned x_base,
    unsigned y_base,
    unsigned distance,
    uint16_t multiplier_log,
    uint64_t byte_count)
{
    // Every pair in a transform range uses the same fixed multiplier.  Keep
    // its eight nibble rows prepared across the independent pairs instead of
    // broadcasting them again for every shard pair.  The AVX2 product still
    // uses the established byte order and arithmetic within each butterfly.
    static const uint16_t kZeroSkew = 65535;
    if (multiplier_log == kZeroSkew)
    {
        for (unsigned i = 0; i < distance; ++i)
            AVX2XorMemory(work[y_base + i], work[x_base + i], byte_count);
        return;
    }
    const FF16Table& table = FF16Tables[multiplier_log];
#if defined(LEO2_GFNI_VARIANT)
    const __m256i low_tables[4] = {
        BroadcastAffine(table.block[0]), BroadcastAffine(table.block[1]),
        BroadcastAffine(table.block[2]), BroadcastAffine(table.block[3])
    };
    // ProductVectors ignores its high_tables operand in this variant.
    const __m256i* const high_tables = low_tables;
#else
    const __m256i low_tables[4] = {
        BroadcastTable(table.low[0]), BroadcastTable(table.low[1]),
        BroadcastTable(table.low[2]), BroadcastTable(table.low[3])
    };
    const __m256i high_tables[4] = {
        BroadcastTable(table.high[0]), BroadcastTable(table.high[1]),
        BroadcastTable(table.high[2]), BroadcastTable(table.high[3])
    };
#endif
    for (unsigned i = 0; i < distance; ++i)
    {
        AVX2FF16Butterfly2Prepared<Inverse>(
            work[x_base + i], work[y_base + i], multiplier_log,
            low_tables, high_tables, byte_count);
    }
}
#endif

template<bool Inverse>
static void AVX2FF16Butterfly4Range(
    void* const* work, unsigned distance,
    uint16_t log01, uint16_t log23, uint16_t log02,
    uint64_t byte_count, bool prefer_fused)
{
#if defined(LEO2_GFNI_VARIANT)
    // Affine multipliers cost four re-broadcast vectors instead of eight held
    // ones, so the fused radix-four group fits at every size and the caller's
    // prefer_fused hint, which exists to bound nibble-table register pressure,
    // no longer applies to this backend.
    (void)prefer_fused;
    for (unsigned i = 0; i < distance; ++i)
    {
        AVX2FF16Butterfly4<Inverse>(
            work[i], work[i + distance],
            work[i + distance * 2U], work[i + distance * 3U],
            log01, log23, log02, byte_count);
    }
#elif defined(LEO2_AVX512_VARIANT)
    for (unsigned i = 0; i < distance; ++i)
    {
        void* const value0 = work[i];
        void* const value1 = work[i + distance];
        void* const value2 = work[i + distance * 2U];
        void* const value3 = work[i + distance * 3U];
        if (prefer_fused)
            AVX2FF16Butterfly4<Inverse>(value0, value1, value2, value3,
                log01, log23, log02, byte_count);
        else
            AVX2FF16Butterfly4Split<Inverse>(
                value0, value1, value2, value3,
                log01, log23, log02, byte_count);
    }
#else
    // The Butterfly4Range contract makes all coordinates disjoint, so the
    // split radix-four layers are independent across i.  Execute each layer
    // as a range so its fixed tables remain prepared, then preserve the
    // inverse/forward layer order for every coordinate.  AVX-512VL retains
    // the fused/per-pair schedule above: its expanded register file has a
    // different measured crossover and does not need this AVX2 workaround.
    if (prefer_fused)
    {
        for (unsigned i = 0; i < distance; ++i)
        {
            AVX2FF16Butterfly4<Inverse>(
                work[i], work[i + distance],
                work[i + distance * 2U], work[i + distance * 3U],
                log01, log23, log02, byte_count);
        }
        return;
    }
    if (Inverse)
    {
        AVX2FF16Butterfly2RangePrepared<true>(
            work, 0, distance, distance, log01, byte_count);
        AVX2FF16Butterfly2RangePrepared<true>(
            work, distance * 2U, distance * 3U,
            distance, log23, byte_count);
        AVX2FF16Butterfly2RangePrepared<true>(
            work, 0, distance * 2U, distance, log02, byte_count);
        AVX2FF16Butterfly2RangePrepared<true>(
            work, distance, distance * 3U,
            distance, log02, byte_count);
    }
    else
    {
        AVX2FF16Butterfly2RangePrepared<false>(
            work, 0, distance * 2U, distance, log02, byte_count);
        AVX2FF16Butterfly2RangePrepared<false>(
            work, distance, distance * 3U,
            distance, log02, byte_count);
        AVX2FF16Butterfly2RangePrepared<false>(
            work, 0, distance, distance, log01, byte_count);
        AVX2FF16Butterfly2RangePrepared<false>(
            work, distance * 2U, distance * 3U,
            distance, log23, byte_count);
    }
#endif
}

static void AVX2FF16IFFTButterfly4Range(
    void* const* work, unsigned distance,
    uint16_t log01, uint16_t log23, uint16_t log02,
    uint64_t byte_count, bool prefer_fused)
{
    AVX2FF16Butterfly4Range<true>(work, distance,
        log01, log23, log02, byte_count, prefer_fused);
}

static void AVX2FF16FFTButterfly4Range(
    void* const* work, unsigned distance,
    uint16_t log01, uint16_t log23, uint16_t log02,
    uint64_t byte_count, bool prefer_fused)
{
    AVX2FF16Butterfly4Range<false>(work, distance,
        log01, log23, log02, byte_count, prefer_fused);
}
#endif // LEO_HAS_FF16

#if defined(LEO_HAS_FF8) && !defined(LEO2_AVX512_VARIANT)
#if defined(_MSC_VER)
# define LEO2_AVX2_DIRECT_NOINLINE __declspec(noinline)
#elif defined(__GNUC__) || defined(__clang__)
# define LEO2_AVX2_DIRECT_NOINLINE __attribute__((noinline))
#else
# define LEO2_AVX2_DIRECT_NOINLINE
#endif

static LEO2_AVX2_DIRECT_NOINLINE void
AVX2FF8MultiplyAddOutputPair(
    void* destination0_pointer,
    void* destination1_pointer,
    const void* source_pointer,
    uint16_t log0,
    uint16_t log1,
    uint64_t byte_count)
{
    uint8_t* const destination0 =
        static_cast<uint8_t*>(destination0_pointer);
    uint8_t* const destination1 =
        static_cast<uint8_t*>(destination1_pointer);
    const uint8_t* const source =
        static_cast<const uint8_t*>(source_pointer);
    const FF8NibbleTable& table0 = FF8Tables[log0];
    const FF8NibbleTable& table1 = FF8Tables[log1];
    const __m256i low0 = BroadcastTable(table0.low);
#if !defined(LEO2_GFNI_VARIANT)
    const __m256i high0 = BroadcastTable(table0.high);
#endif
    const __m256i low1 = BroadcastTable(table1.low);
#if !defined(LEO2_GFNI_VARIANT)
    const __m256i high1 = BroadcastTable(table1.high);
    const __m256i nibble_mask = _mm256_set1_epi8(0x0f);
#endif
    uint64_t offset = 0;
    while (byte_count - offset >= 32)
    {
        const __m256i data = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(source + offset));
#if defined(LEO2_GFNI_VARIANT)
        __m256i product = _mm256_gf2p8affine_epi64_epi8(data, low0, 0);
#else
        const __m256i low_indices = _mm256_and_si256(data, nibble_mask);
        const __m256i high_indices = _mm256_and_si256(
            _mm256_srli_epi64(data, 4), nibble_mask);
        __m256i product = _mm256_xor_si256(
            _mm256_shuffle_epi8(low0, low_indices),
            _mm256_shuffle_epi8(high0, high_indices));
#endif
        product = _mm256_xor_si256(product, _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(destination0 + offset)));
        _mm256_storeu_si256(reinterpret_cast<__m256i*>(
            destination0 + offset), product);
#if defined(LEO2_GFNI_VARIANT)
        product = _mm256_gf2p8affine_epi64_epi8(data, low1, 0);
#else
        product = _mm256_xor_si256(
            _mm256_shuffle_epi8(low1, low_indices),
            _mm256_shuffle_epi8(high1, high_indices));
#endif
        product = _mm256_xor_si256(product, _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(destination1 + offset)));
        _mm256_storeu_si256(reinterpret_cast<__m256i*>(
            destination1 + offset), product);
        offset += 32;
    }
    while (offset < byte_count)
    {
        const uint8_t value = source[offset];
        destination0[offset] ^= FF8Product(log0, value);
        destination1[offset] ^= FF8Product(log1, value);
        ++offset;
    }
}

static LEO2_AVX2_DIRECT_NOINLINE void
AVX2FF8MultiplyAddOutputGroup2(
    void* const* destination_pointers,
    const void* source_pointer,
    const uint16_t* multiplier_logs,
    uint64_t byte_count)
{
    const uint16_t log0 = multiplier_logs[0];
    const uint16_t log1 = multiplier_logs[1];
    if (log0 != UINT16_MAX && log1 != UINT16_MAX)
    {
        AVX2FF8MultiplyAddOutputPair(
            destination_pointers[0], destination_pointers[1],
            source_pointer, log0, log1, byte_count);
        return;
    }
    if (log0 != UINT16_MAX)
        AVX2FF8MultiplyAdd(destination_pointers[0],
            source_pointer, log0, byte_count);
    if (log1 != UINT16_MAX)
        AVX2FF8MultiplyAdd(destination_pointers[1],
            source_pointer, log1, byte_count);
}

static LEO2_AVX2_DIRECT_NOINLINE void
AVX2FF8MultiplyAddOutputGroup4(
    void* const* destination_pointers,
    const void* source_pointer,
    const uint16_t* multiplier_logs,
    uint64_t byte_count)
{
    /*
        Direct-repair source rows are normally dense: every surviving source
        contributes to every missing original.  Keep those eight nibble
        tables in explicitly named registers.  GCC otherwise materializes
        the indexed table arrays on the stack before reloading them into the
        sixteen-register AVX2 file for each source row.

        A selected parity row may omit the output it initialized.  Such sparse
        rows are uncommon and use the mature single-output kernel rather than
        reserving eight vector slots for inactive tables.
    */
    const bool all_active =
        multiplier_logs[0] != UINT16_MAX &&
        multiplier_logs[1] != UINT16_MAX &&
        multiplier_logs[2] != UINT16_MAX &&
        multiplier_logs[3] != UINT16_MAX;
    if (!all_active)
    {
        void* paired_destination0 = NULL;
        void* paired_destination1 = NULL;
        uint16_t paired_log0 = UINT16_MAX;
        uint16_t paired_log1 = UINT16_MAX;
        unsigned paired_count = 0;
        for (unsigned output = 0; output < 4; ++output)
        {
            const uint16_t log = multiplier_logs[output];
            if (log != UINT16_MAX)
            {
                if (paired_count < 2)
                {
                    if (paired_count++ == 0)
                    {
                        paired_destination0 = destination_pointers[output];
                        paired_log0 = log;
                    }
                    else
                    {
                        paired_destination1 = destination_pointers[output];
                        paired_log1 = log;
                    }
                }
                else
                {
                    AVX2FF8MultiplyAdd(
                        destination_pointers[output], source_pointer,
                        log, byte_count);
                }
            }
        }
        if (paired_count == 2)
        {
            AVX2FF8MultiplyAddOutputPair(
                paired_destination0, paired_destination1,
                source_pointer, paired_log0, paired_log1, byte_count);
        }
        else if (paired_count == 1)
        {
            AVX2FF8MultiplyAdd(paired_destination0,
                source_pointer, paired_log0, byte_count);
        }
        return;
    }

    uint8_t* const destination0 =
        static_cast<uint8_t*>(destination_pointers[0]);
    uint8_t* const destination1 =
        static_cast<uint8_t*>(destination_pointers[1]);
    uint8_t* const destination2 =
        static_cast<uint8_t*>(destination_pointers[2]);
    uint8_t* const destination3 =
        static_cast<uint8_t*>(destination_pointers[3]);
    const uint8_t* const source =
        static_cast<const uint8_t*>(source_pointer);
    const uint16_t log0 = multiplier_logs[0];
    const uint16_t log1 = multiplier_logs[1];
    const uint16_t log2 = multiplier_logs[2];
    const uint16_t log3 = multiplier_logs[3];
    const FF8NibbleTable& table0 = FF8Tables[log0];
    const FF8NibbleTable& table1 = FF8Tables[log1];
    const FF8NibbleTable& table2 = FF8Tables[log2];
    const FF8NibbleTable& table3 = FF8Tables[log3];
    const __m256i low0 = BroadcastTable(table0.low);
#if !defined(LEO2_GFNI_VARIANT)
    const __m256i high0 = BroadcastTable(table0.high);
#endif
    const __m256i low1 = BroadcastTable(table1.low);
#if !defined(LEO2_GFNI_VARIANT)
    const __m256i high1 = BroadcastTable(table1.high);
#endif
    const __m256i low2 = BroadcastTable(table2.low);
#if !defined(LEO2_GFNI_VARIANT)
    const __m256i high2 = BroadcastTable(table2.high);
#endif
    const __m256i low3 = BroadcastTable(table3.low);
#if !defined(LEO2_GFNI_VARIANT)
    const __m256i high3 = BroadcastTable(table3.high);
    const __m256i nibble_mask = _mm256_set1_epi8(0x0f);
#endif
    uint64_t offset = 0;
    while (byte_count - offset >= 32)
    {
        const __m256i data = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(source + offset));
#if !defined(LEO2_GFNI_VARIANT)
        const __m256i low_indices = _mm256_and_si256(data, nibble_mask);
        const __m256i high_indices = _mm256_and_si256(
            _mm256_srli_epi64(data, 4), nibble_mask);
#endif
#if defined(LEO2_GFNI_VARIANT)
        __m256i product = _mm256_gf2p8affine_epi64_epi8(data, low0, 0);
#else
        __m256i product = _mm256_xor_si256(
            _mm256_shuffle_epi8(low0, low_indices),
            _mm256_shuffle_epi8(high0, high_indices));
#endif
        product = _mm256_xor_si256(product, _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(destination0 + offset)));
        _mm256_storeu_si256(reinterpret_cast<__m256i*>(
            destination0 + offset), product);
#if defined(LEO2_GFNI_VARIANT)
        product = _mm256_gf2p8affine_epi64_epi8(data, low1, 0);
#else
        product = _mm256_xor_si256(
            _mm256_shuffle_epi8(low1, low_indices),
            _mm256_shuffle_epi8(high1, high_indices));
#endif
        product = _mm256_xor_si256(product, _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(destination1 + offset)));
        _mm256_storeu_si256(reinterpret_cast<__m256i*>(
            destination1 + offset), product);
#if defined(LEO2_GFNI_VARIANT)
        product = _mm256_gf2p8affine_epi64_epi8(data, low2, 0);
#else
        product = _mm256_xor_si256(
            _mm256_shuffle_epi8(low2, low_indices),
            _mm256_shuffle_epi8(high2, high_indices));
#endif
        product = _mm256_xor_si256(product, _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(destination2 + offset)));
        _mm256_storeu_si256(reinterpret_cast<__m256i*>(
            destination2 + offset), product);
#if defined(LEO2_GFNI_VARIANT)
        product = _mm256_gf2p8affine_epi64_epi8(data, low3, 0);
#else
        product = _mm256_xor_si256(
            _mm256_shuffle_epi8(low3, low_indices),
            _mm256_shuffle_epi8(high3, high_indices));
#endif
        product = _mm256_xor_si256(product, _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(destination3 + offset)));
        _mm256_storeu_si256(reinterpret_cast<__m256i*>(
            destination3 + offset), product);
        offset += 32;
    }
    while (offset < byte_count)
    {
        const uint8_t value = source[offset];
        destination0[offset] ^= FF8Product(log0, value);
        destination1[offset] ^= FF8Product(log1, value);
        destination2[offset] ^= FF8Product(log2, value);
        destination3[offset] ^= FF8Product(log3, value);
        ++offset;
    }
}

#if !defined(LEO2_GFNI_VARIANT)
#if defined(_MSC_VER)
#define LEO2_AVX2_GROUP35_NOINLINE __declspec(noinline)
#elif (defined(__GNUC__) || defined(__clang__)) && defined(__ELF__)
#define LEO2_AVX2_GROUP35_NOINLINE \
    __attribute__((noinline, section(".text.leo2_avx2_direct_group35")))
#elif defined(__GNUC__) || defined(__clang__)
#define LEO2_AVX2_GROUP35_NOINLINE __attribute__((noinline))
#else
#define LEO2_AVX2_GROUP35_NOINLINE
#endif

static LEO2_AVX2_GROUP35_NOINLINE void
AVX2FF8MultiplyAddOutputGroup3(
    void* const* destination_pointers,
    const void* source_pointer,
    const uint16_t* multiplier_logs,
    uint64_t byte_count)
{
    const bool all_active =
        multiplier_logs[0] != UINT16_MAX &&
        multiplier_logs[1] != UINT16_MAX &&
        multiplier_logs[2] != UINT16_MAX;
    if (!all_active)
    {
        AVX2FF8MultiplyAddOutputGroup2(
            destination_pointers, source_pointer,
            multiplier_logs, byte_count);
        if (multiplier_logs[2] != UINT16_MAX)
        {
            AVX2FF8MultiplyAdd(
                destination_pointers[2], source_pointer,
                multiplier_logs[2], byte_count);
        }
        return;
    }

    uint8_t* const destination0 =
        static_cast<uint8_t*>(destination_pointers[0]);
    uint8_t* const destination1 =
        static_cast<uint8_t*>(destination_pointers[1]);
    uint8_t* const destination2 =
        static_cast<uint8_t*>(destination_pointers[2]);
    const uint8_t* const source =
        static_cast<const uint8_t*>(source_pointer);
    const uint16_t log0 = multiplier_logs[0];
    const uint16_t log1 = multiplier_logs[1];
    const uint16_t log2 = multiplier_logs[2];
    const FF8NibbleTable& table0 = FF8Tables[log0];
    const FF8NibbleTable& table1 = FF8Tables[log1];
    const FF8NibbleTable& table2 = FF8Tables[log2];
    const __m256i low0 = BroadcastTable(table0.low);
    const __m256i high0 = BroadcastTable(table0.high);
    const __m256i low1 = BroadcastTable(table1.low);
    const __m256i high1 = BroadcastTable(table1.high);
    const __m256i low2 = BroadcastTable(table2.low);
    const __m256i high2 = BroadcastTable(table2.high);
    const __m256i nibble_mask = _mm256_set1_epi8(0x0f);
    uint64_t offset = 0;
    while (byte_count - offset >= 32)
    {
        const __m256i data = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(source + offset));
        const __m256i low_indices =
            _mm256_and_si256(data, nibble_mask);
        const __m256i high_indices = _mm256_and_si256(
            _mm256_srli_epi64(data, 4), nibble_mask);
        __m256i product = _mm256_xor_si256(
            _mm256_shuffle_epi8(low0, low_indices),
            _mm256_shuffle_epi8(high0, high_indices));
        product = _mm256_xor_si256(product, _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(destination0 + offset)));
        _mm256_storeu_si256(reinterpret_cast<__m256i*>(
            destination0 + offset), product);
        product = _mm256_xor_si256(
            _mm256_shuffle_epi8(low1, low_indices),
            _mm256_shuffle_epi8(high1, high_indices));
        product = _mm256_xor_si256(product, _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(destination1 + offset)));
        _mm256_storeu_si256(reinterpret_cast<__m256i*>(
            destination1 + offset), product);
        product = _mm256_xor_si256(
            _mm256_shuffle_epi8(low2, low_indices),
            _mm256_shuffle_epi8(high2, high_indices));
        product = _mm256_xor_si256(product, _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(destination2 + offset)));
        _mm256_storeu_si256(reinterpret_cast<__m256i*>(
            destination2 + offset), product);
        offset += 32;
    }
    while (offset < byte_count)
    {
        const uint8_t value = source[offset];
        destination0[offset] ^= FF8Product(log0, value);
        destination1[offset] ^= FF8Product(log1, value);
        destination2[offset] ^= FF8Product(log2, value);
        ++offset;
    }
}

static LEO2_AVX2_GROUP35_NOINLINE void
AVX2FF8MultiplyAddOutputGroup5(
    void* const* destination_pointers,
    const void* source_pointer,
    const uint16_t* multiplier_logs,
    uint64_t byte_count)
{
    const bool all_active =
        multiplier_logs[0] != UINT16_MAX &&
        multiplier_logs[1] != UINT16_MAX &&
        multiplier_logs[2] != UINT16_MAX &&
        multiplier_logs[3] != UINT16_MAX &&
        multiplier_logs[4] != UINT16_MAX;
    if (!all_active)
    {
        AVX2FF8MultiplyAddOutputGroup4(
            destination_pointers, source_pointer,
            multiplier_logs, byte_count);
        if (multiplier_logs[4] != UINT16_MAX)
        {
            AVX2FF8MultiplyAdd(
                destination_pointers[4], source_pointer,
                multiplier_logs[4], byte_count);
        }
        return;
    }

    uint8_t* const destination0 =
        static_cast<uint8_t*>(destination_pointers[0]);
    uint8_t* const destination1 =
        static_cast<uint8_t*>(destination_pointers[1]);
    uint8_t* const destination2 =
        static_cast<uint8_t*>(destination_pointers[2]);
    uint8_t* const destination3 =
        static_cast<uint8_t*>(destination_pointers[3]);
    uint8_t* const destination4 =
        static_cast<uint8_t*>(destination_pointers[4]);
    const uint8_t* const source =
        static_cast<const uint8_t*>(source_pointer);
    const uint16_t log0 = multiplier_logs[0];
    const uint16_t log1 = multiplier_logs[1];
    const uint16_t log2 = multiplier_logs[2];
    const uint16_t log3 = multiplier_logs[3];
    const uint16_t log4 = multiplier_logs[4];
    const FF8NibbleTable& table0 = FF8Tables[log0];
    const FF8NibbleTable& table1 = FF8Tables[log1];
    const FF8NibbleTable& table2 = FF8Tables[log2];
    const FF8NibbleTable& table3 = FF8Tables[log3];
    const FF8NibbleTable& table4 = FF8Tables[log4];
    const __m256i low0 = BroadcastTable(table0.low);
    const __m256i high0 = BroadcastTable(table0.high);
    const __m256i low1 = BroadcastTable(table1.low);
    const __m256i high1 = BroadcastTable(table1.high);
    const __m256i low2 = BroadcastTable(table2.low);
    const __m256i high2 = BroadcastTable(table2.high);
    const __m256i low3 = BroadcastTable(table3.low);
    const __m256i high3 = BroadcastTable(table3.high);
    const __m256i low4 = BroadcastTable(table4.low);
    const __m256i high4 = BroadcastTable(table4.high);
    const __m256i nibble_mask = _mm256_set1_epi8(0x0f);
    uint64_t offset = 0;
    while (byte_count - offset >= 32)
    {
        const __m256i data = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(source + offset));
        const __m256i low_indices =
            _mm256_and_si256(data, nibble_mask);
        const __m256i high_indices = _mm256_and_si256(
            _mm256_srli_epi64(data, 4), nibble_mask);
        __m256i product = _mm256_xor_si256(
            _mm256_shuffle_epi8(low0, low_indices),
            _mm256_shuffle_epi8(high0, high_indices));
        product = _mm256_xor_si256(product, _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(destination0 + offset)));
        _mm256_storeu_si256(reinterpret_cast<__m256i*>(
            destination0 + offset), product);
        product = _mm256_xor_si256(
            _mm256_shuffle_epi8(low1, low_indices),
            _mm256_shuffle_epi8(high1, high_indices));
        product = _mm256_xor_si256(product, _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(destination1 + offset)));
        _mm256_storeu_si256(reinterpret_cast<__m256i*>(
            destination1 + offset), product);
        product = _mm256_xor_si256(
            _mm256_shuffle_epi8(low2, low_indices),
            _mm256_shuffle_epi8(high2, high_indices));
        product = _mm256_xor_si256(product, _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(destination2 + offset)));
        _mm256_storeu_si256(reinterpret_cast<__m256i*>(
            destination2 + offset), product);
        product = _mm256_xor_si256(
            _mm256_shuffle_epi8(low3, low_indices),
            _mm256_shuffle_epi8(high3, high_indices));
        product = _mm256_xor_si256(product, _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(destination3 + offset)));
        _mm256_storeu_si256(reinterpret_cast<__m256i*>(
            destination3 + offset), product);
        product = _mm256_xor_si256(
            _mm256_shuffle_epi8(low4, low_indices),
            _mm256_shuffle_epi8(high4, high_indices));
        product = _mm256_xor_si256(product, _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(destination4 + offset)));
        _mm256_storeu_si256(reinterpret_cast<__m256i*>(
            destination4 + offset), product);
        offset += 32;
    }
    while (offset < byte_count)
    {
        const uint8_t value = source[offset];
        destination0[offset] ^= FF8Product(log0, value);
        destination1[offset] ^= FF8Product(log1, value);
        destination2[offset] ^= FF8Product(log2, value);
        destination3[offset] ^= FF8Product(log3, value);
        destination4[offset] ^= FF8Product(log4, value);
        ++offset;
    }
}

#undef LEO2_AVX2_GROUP35_NOINLINE

#if (defined(__GNUC__) || defined(__clang__)) && !defined(_MSC_VER)
#if defined(__ELF__)
#define LEO2_AVX2_GROUP6_NOINLINE \
    __attribute__((noinline, section(".text.leo2_avx2_direct_group6")))
#else
#define LEO2_AVX2_GROUP6_NOINLINE __attribute__((noinline))
#endif

static const uint8_t AVX2FF8NibbleMaskByte = 0x0f;

static LEO2_AVX2_GROUP6_NOINLINE void
AVX2FF8MultiplyAddOutputGroup6(
    void* const* destination_pointers,
    const void* source_pointer,
    const uint16_t* multiplier_logs,
    uint64_t byte_count)
{
    const bool all_active =
        multiplier_logs[0] != UINT16_MAX &&
        multiplier_logs[1] != UINT16_MAX &&
        multiplier_logs[2] != UINT16_MAX &&
        multiplier_logs[3] != UINT16_MAX &&
        multiplier_logs[4] != UINT16_MAX &&
        multiplier_logs[5] != UINT16_MAX;
    if (!all_active)
    {
        AVX2FF8MultiplyAddOutputGroup4(
            destination_pointers, source_pointer,
            multiplier_logs, byte_count);
        AVX2FF8MultiplyAddOutputGroup2(
            destination_pointers + 4, source_pointer,
            multiplier_logs + 4, byte_count);
        return;
    }

    uint8_t* const destination0 =
        static_cast<uint8_t*>(destination_pointers[0]);
    uint8_t* const destination1 =
        static_cast<uint8_t*>(destination_pointers[1]);
    uint8_t* const destination2 =
        static_cast<uint8_t*>(destination_pointers[2]);
    uint8_t* const destination3 =
        static_cast<uint8_t*>(destination_pointers[3]);
    uint8_t* const destination4 =
        static_cast<uint8_t*>(destination_pointers[4]);
    uint8_t* const destination5 =
        static_cast<uint8_t*>(destination_pointers[5]);
    const uint8_t* const source =
        static_cast<const uint8_t*>(source_pointer);
    const uint16_t log0 = multiplier_logs[0];
    const uint16_t log1 = multiplier_logs[1];
    const uint16_t log2 = multiplier_logs[2];
    const uint16_t log3 = multiplier_logs[3];
    const uint16_t log4 = multiplier_logs[4];
    const uint16_t log5 = multiplier_logs[5];
    const FF8NibbleTable& table0 = FF8Tables[log0];
    const FF8NibbleTable& table1 = FF8Tables[log1];
    const FF8NibbleTable& table2 = FF8Tables[log2];
    const FF8NibbleTable& table3 = FF8Tables[log3];
    const FF8NibbleTable& table4 = FF8Tables[log4];
    const FF8NibbleTable& table5 = FF8Tables[log5];
    const __m256i low0 = BroadcastTable(table0.low);
    const __m256i high0 = BroadcastTable(table0.high);
    const __m256i low1 = BroadcastTable(table1.low);
    const __m256i high1 = BroadcastTable(table1.high);
    const __m256i low2 = BroadcastTable(table2.low);
    const __m256i high2 = BroadcastTable(table2.high);
    const __m256i low3 = BroadcastTable(table3.low);
    const __m256i high3 = BroadcastTable(table3.high);
    const __m256i low4 = BroadcastTable(table4.low);
    const __m256i high4 = BroadcastTable(table4.high);
    const __m256i low5 = BroadcastTable(table5.low);
    const __m256i high5 = BroadcastTable(table5.high);
    // The dispatcher admits only complete AVX2 vectors to this dense kernel.
    uint64_t offset = 0;
    while (offset < byte_count)
    {
        __m256i high_indices = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(source + offset));
        __m256i low_indices;
        __m256i generated_mask;
        __asm__(
            "vpbroadcastb %3, %2\n\t"
            "vpand %2, %1, %0\n\t"
            "vpsrlq $4, %1, %1\n\t"
            "vpand %2, %1, %1"
            : "=&x"(low_indices), "+x"(high_indices),
              "=&x"(generated_mask)
            : "m"(AVX2FF8NibbleMaskByte));
        __m256i product = _mm256_xor_si256(
            _mm256_shuffle_epi8(low0, low_indices),
            _mm256_shuffle_epi8(high0, high_indices));
        product = _mm256_xor_si256(product, _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(destination0 + offset)));
        _mm256_storeu_si256(reinterpret_cast<__m256i*>(
            destination0 + offset), product);
        product = _mm256_xor_si256(
            _mm256_shuffle_epi8(low1, low_indices),
            _mm256_shuffle_epi8(high1, high_indices));
        product = _mm256_xor_si256(product, _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(destination1 + offset)));
        _mm256_storeu_si256(reinterpret_cast<__m256i*>(
            destination1 + offset), product);
        product = _mm256_xor_si256(
            _mm256_shuffle_epi8(low2, low_indices),
            _mm256_shuffle_epi8(high2, high_indices));
        product = _mm256_xor_si256(product, _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(destination2 + offset)));
        _mm256_storeu_si256(reinterpret_cast<__m256i*>(
            destination2 + offset), product);
        product = _mm256_xor_si256(
            _mm256_shuffle_epi8(low3, low_indices),
            _mm256_shuffle_epi8(high3, high_indices));
        product = _mm256_xor_si256(product, _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(destination3 + offset)));
        _mm256_storeu_si256(reinterpret_cast<__m256i*>(
            destination3 + offset), product);
        product = _mm256_xor_si256(
            _mm256_shuffle_epi8(low4, low_indices),
            _mm256_shuffle_epi8(high4, high_indices));
        product = _mm256_xor_si256(product, _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(destination4 + offset)));
        _mm256_storeu_si256(reinterpret_cast<__m256i*>(
            destination4 + offset), product);
        product = _mm256_xor_si256(
            _mm256_shuffle_epi8(low5, low_indices),
            _mm256_shuffle_epi8(high5, high_indices));
        product = _mm256_xor_si256(product, _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(destination5 + offset)));
        _mm256_storeu_si256(reinterpret_cast<__m256i*>(
            destination5 + offset), product);
        offset += 32;
    }
}

#undef LEO2_AVX2_GROUP6_NOINLINE
#endif
#endif

#undef LEO2_AVX2_DIRECT_NOINLINE

static void AVX2FF8MultiplyAddOutputs(
    void* const* destinations,
    const void* source,
    const uint16_t* multiplier_logs,
    uint32_t output_count,
    uint64_t byte_count)
{
    if (output_count == 2)
    {
        AVX2FF8MultiplyAddOutputGroup2(
            destinations, source, multiplier_logs, byte_count);
        return;
    }
    if (output_count == 3)
    {
#if !defined(LEO2_GFNI_VARIANT)
        AVX2FF8MultiplyAddOutputGroup3(
            destinations, source, multiplier_logs, byte_count);
#else
        AVX2FF8MultiplyAddOutputGroup2(
            destinations, source, multiplier_logs, byte_count);
        if (multiplier_logs[2] != UINT16_MAX)
        {
            AVX2FF8MultiplyAdd(destinations[2], source,
                multiplier_logs[2], byte_count);
        }
#endif
        return;
    }
    if (output_count == 4)
    {
        AVX2FF8MultiplyAddOutputGroup4(
            destinations, source, multiplier_logs, byte_count);
        return;
    }
    if (output_count < 5 || output_count > 8)
        return;

#if !defined(LEO2_GFNI_VARIANT)
    if (output_count == 5)
    {
        AVX2FF8MultiplyAddOutputGroup5(
            destinations, source, multiplier_logs, byte_count);
        return;
    }
    if (output_count == 7)
    {
        AVX2FF8MultiplyAddOutputGroup4(
            destinations, source, multiplier_logs, byte_count);
        AVX2FF8MultiplyAddOutputGroup3(
            destinations + 4, source, multiplier_logs + 4, byte_count);
        return;
    }
#if (defined(__GNUC__) || defined(__clang__)) && !defined(_MSC_VER)
    if (output_count == 6 &&
        byte_count >= UINT64_C(2048) &&
        (byte_count & UINT64_C(31)) == 0)
    {
        AVX2FF8MultiplyAddOutputGroup6(
            destinations, source, multiplier_logs, byte_count);
        return;
    }
#endif
#endif

    // Compose the established zero-spill four- and two-output kernels and
    // handle an optional final output.  This rereads each source at most
    // three times, versus once per missing original in the output-major
    // executor.  Eight outputs retain two regular four-output passes because
    // sixteen nibble-table vectors do not fit alongside the source indices
    // and destinations in AVX2's sixteen-register file.
    AVX2FF8MultiplyAddOutputGroup4(
        destinations, source, multiplier_logs, byte_count);
    const uint32_t remaining = output_count - 4;
    if (remaining == 4)
    {
        AVX2FF8MultiplyAddOutputGroup4(
            destinations + 4, source, multiplier_logs + 4, byte_count);
        return;
    }
    if (remaining >= 2)
    {
        AVX2FF8MultiplyAddOutputGroup2(
            destinations + 4, source, multiplier_logs + 4, byte_count);
    }
    if ((remaining & 1U) != 0 &&
        multiplier_logs[output_count - 1] != UINT16_MAX)
    {
        AVX2FF8MultiplyAdd(destinations[output_count - 1], source,
            multiplier_logs[output_count - 1], byte_count);
    }
}

#if !defined(LEO2_AVX512_VARIANT) && !defined(LEO2_GFNI_VARIANT)
#if defined(_MSC_VER)
# define LEO2_AVX2_TINY_GROUP_NOINLINE __declspec(noinline)
#elif defined(__GNUC__) || defined(__clang__)
# define LEO2_AVX2_TINY_GROUP_NOINLINE __attribute__((noinline))
#else
# define LEO2_AVX2_TINY_GROUP_NOINLINE
#endif

template<uint32_t Width>
static LEO_FORCE_INLINE __m128i
AVX2FF8TinyLoadXMM(const uint8_t* source)
{
    static_assert(
        Width == 1 || Width == 2 || Width == 4 ||
            Width == 8 || Width == 16,
        "tiny XMM load width must be 1, 2, 4, 8, or 16 bytes");
    if (Width == 16)
    {
        return _mm_loadu_si128(
            reinterpret_cast<const __m128i*>(source));
    }
    if (Width == 8)
    {
        return _mm_loadl_epi64(
            reinterpret_cast<const __m128i*>(source));
    }
    uint32_t value = 0;
    memcpy(&value, source, Width);
    return _mm_cvtsi32_si128(static_cast<int>(value));
}

template<uint32_t Width>
static LEO_FORCE_INLINE void AVX2FF8TinyStoreXMM(
    uint8_t* destination,
    __m128i value)
{
    static_assert(
        Width == 1 || Width == 2 || Width == 4 ||
            Width == 8 || Width == 16,
        "tiny XMM store width must be 1, 2, 4, 8, or 16 bytes");
    if (Width == 16)
    {
        _mm_storeu_si128(
            reinterpret_cast<__m128i*>(destination), value);
        return;
    }
    if (Width == 8)
    {
        _mm_storel_epi64(
            reinterpret_cast<__m128i*>(destination), value);
        return;
    }
    const uint32_t word =
        static_cast<uint32_t>(_mm_cvtsi128_si32(value));
    memcpy(destination, &word, Width);
}

template<
    uint32_t OutputCount,
    uint32_t Width,
    uint32_t StaticSourceCount = 0>
static LEO_FORCE_INLINE void
AVX2FF8EncodeOutputGroupTinyXMM(
    const uint8_t* const* sources,
    uint32_t source_count,
    uint8_t* const* destinations,
    const uint8_t* const* rows,
    uint64_t offset)
{
    const uint32_t source_limit =
        StaticSourceCount != 0 ? StaticSourceCount : source_count;
    __m128i accumulators[OutputCount];
    for (uint32_t output = 0; output < OutputCount; ++output)
        accumulators[output] = _mm_setzero_si128();
    const __m128i nibble_mask = _mm_set1_epi8(0x0f);

    for (uint32_t source = 0; source < source_limit; ++source)
    {
        const __m128i data =
            AVX2FF8TinyLoadXMM<Width>(sources[source] + offset);
        const __m128i low_indices =
            _mm_and_si128(data, nibble_mask);
        const __m128i high_indices = _mm_and_si128(
            _mm_srli_epi64(data, 4), nibble_mask);
        for (uint32_t output = 0; output < OutputCount; ++output)
        {
            const uint16_t log = rows[output][source];
            __m128i product = data;
            if (log != 0)
            {
                const FF8NibbleTable& table = FF8Tables[log];
                product = _mm_xor_si128(
                    _mm_shuffle_epi8(
                        _mm_loadu_si128(
                            reinterpret_cast<const __m128i*>(table.low)),
                        low_indices),
                    _mm_shuffle_epi8(
                        _mm_loadu_si128(
                            reinterpret_cast<const __m128i*>(table.high)),
                        high_indices));
            }
            accumulators[output] = _mm_xor_si128(
                accumulators[output], product);
        }
    }
    for (uint32_t output = 0; output < OutputCount; ++output)
    {
        AVX2FF8TinyStoreXMM<Width>(
            destinations[output] + offset,
            accumulators[output]);
    }
}

template<uint32_t OutputCount, uint32_t StaticSourceCount = 0>
static LEO2_AVX2_TINY_GROUP_NOINLINE void
AVX2FF8EncodeOutputGroupTinyImpl(
    const void* const* source_pointers,
    uint32_t source_count,
    void* const* destination_pointers,
    const uint8_t* coefficient_logs,
    uint64_t byte_count)
{
    static_assert(OutputCount >= 1 && OutputCount <= 8,
        "tiny encode output group must contain one through eight outputs");
    const uint8_t* sources[16];
    uint8_t* destinations[OutputCount];
    const uint8_t* rows[OutputCount];
    const uint32_t source_limit =
        StaticSourceCount != 0 ? StaticSourceCount : source_count;
    for (uint32_t source = 0; source < source_limit; ++source)
        sources[source] =
            static_cast<const uint8_t*>(source_pointers[source]);
    for (uint32_t output = 0; output < OutputCount; ++output)
    {
        destinations[output] =
            static_cast<uint8_t*>(destination_pointers[output]);
        rows[output] = coefficient_logs +
            static_cast<size_t>(output) * source_limit;
    }

    uint64_t offset = 0;
    while (OutputCount <= 4 && byte_count - offset >= 64)
    {
        __m256i low_accumulators[OutputCount];
        __m256i high_accumulators[OutputCount];
        for (uint32_t output = 0; output < OutputCount; ++output)
        {
            low_accumulators[output] = _mm256_setzero_si256();
            high_accumulators[output] = _mm256_setzero_si256();
        }

        for (uint32_t source = 0; source < source_limit; ++source)
        {
            const __m256i low_data = _mm256_loadu_si256(
                reinterpret_cast<const __m256i*>(
                    sources[source] + offset));
            const __m256i high_data = _mm256_loadu_si256(
                reinterpret_cast<const __m256i*>(
                    sources[source] + offset + 32));
            const __m256i nibble_mask = _mm256_set1_epi8(0x0f);
            const __m256i low_low_indices =
                _mm256_and_si256(low_data, nibble_mask);
            const __m256i low_high_indices = _mm256_and_si256(
                _mm256_srli_epi64(low_data, 4), nibble_mask);
            const __m256i high_low_indices =
                _mm256_and_si256(high_data, nibble_mask);
            const __m256i high_high_indices = _mm256_and_si256(
                _mm256_srli_epi64(high_data, 4), nibble_mask);
            for (uint32_t output = 0; output < OutputCount; ++output)
            {
                const uint16_t log = rows[output][source];
                if (log == 0)
                {
                    low_accumulators[output] = _mm256_xor_si256(
                        low_accumulators[output], low_data);
                    high_accumulators[output] = _mm256_xor_si256(
                        high_accumulators[output], high_data);
                }
                else
                {
                    const FF8NibbleTable& table = FF8Tables[log];
                    const __m256i low_table =
                        BroadcastTable(table.low);
                    const __m256i high_table =
                        BroadcastTable(table.high);
                    low_accumulators[output] = _mm256_xor_si256(
                        low_accumulators[output],
                        _mm256_xor_si256(
                            _mm256_shuffle_epi8(
                                low_table, low_low_indices),
                            _mm256_shuffle_epi8(
                                high_table, low_high_indices)));
                    high_accumulators[output] = _mm256_xor_si256(
                        high_accumulators[output],
                        _mm256_xor_si256(
                            _mm256_shuffle_epi8(
                                low_table, high_low_indices),
                            _mm256_shuffle_epi8(
                                high_table, high_high_indices)));
                }
            }
        }
        for (uint32_t output = 0; output < OutputCount; ++output)
        {
            _mm256_storeu_si256(
                reinterpret_cast<__m256i*>(
                    destinations[output] + offset),
                low_accumulators[output]);
            _mm256_storeu_si256(
                reinterpret_cast<__m256i*>(
                    destinations[output] + offset + 32),
                high_accumulators[output]);
        }
        offset += 64;
    }

    while (byte_count - offset >= 32)
    {
        __m256i accumulators[OutputCount];
        for (uint32_t output = 0; output < OutputCount; ++output)
            accumulators[output] = _mm256_setzero_si256();

        const __m256i nibble_mask = _mm256_set1_epi8(0x0f);
        for (uint32_t source = 0; source < source_limit; ++source)
        {
            const __m256i data = _mm256_loadu_si256(
                reinterpret_cast<const __m256i*>(
                    sources[source] + offset));
            const __m256i low_indices =
                _mm256_and_si256(data, nibble_mask);
            const __m256i high_indices = _mm256_and_si256(
                _mm256_srli_epi64(data, 4), nibble_mask);
            for (uint32_t output = 0; output < OutputCount; ++output)
            {
                const uint16_t log = rows[output][source];
                if (log == 0)
                {
                    accumulators[output] = _mm256_xor_si256(
                        accumulators[output], data);
                }
                else
                {
                    const FF8NibbleTable& table = FF8Tables[log];
                    const __m256i low_table =
                        BroadcastTable(table.low);
                    const __m256i high_table =
                        BroadcastTable(table.high);
                    accumulators[output] = _mm256_xor_si256(
                        accumulators[output],
                        _mm256_xor_si256(
                            _mm256_shuffle_epi8(
                                low_table, low_indices),
                            _mm256_shuffle_epi8(
                                high_table, high_indices)));
                }
            }
        }
        for (uint32_t output = 0; output < OutputCount; ++output)
        {
            _mm256_storeu_si256(
                reinterpret_cast<__m256i*>(
                    destinations[output] + offset),
                accumulators[output]);
        }
        offset += 32;
    }

    if (byte_count - offset >= 16)
    {
        AVX2FF8EncodeOutputGroupTinyXMM<
            OutputCount, 16, StaticSourceCount>(
                sources, source_limit, destinations, rows, offset);
        offset += 16;
    }
    if (byte_count - offset >= 8)
    {
        AVX2FF8EncodeOutputGroupTinyXMM<
            OutputCount, 8, StaticSourceCount>(
                sources, source_limit, destinations, rows, offset);
        offset += 8;
    }
    if (byte_count - offset >= 4)
    {
        AVX2FF8EncodeOutputGroupTinyXMM<
            OutputCount, 4, StaticSourceCount>(
                sources, source_limit, destinations, rows, offset);
        offset += 4;
    }

    switch (byte_count - offset)
    {
    case 1:
        AVX2FF8EncodeOutputGroupTinyXMM<
            OutputCount, 1, StaticSourceCount>(
                sources, source_limit, destinations, rows, offset);
        break;
    case 2:
        AVX2FF8EncodeOutputGroupTinyXMM<
            OutputCount, 2, StaticSourceCount>(
                sources, source_limit, destinations, rows, offset);
        break;
    case 3:
        for (uint32_t tail = 0; tail < 3; ++tail)
        {
            for (uint32_t output = 0; output < OutputCount; ++output)
            {
                uint8_t value = 0;
                for (uint32_t source = 0;
                     source < source_limit;
                     ++source)
                {
                    const uint16_t log = rows[output][source];
                    value ^= log == 0
                        ? sources[source][offset + tail]
                        : FF8Product(
                            log, sources[source][offset + tail]);
                }
                destinations[output][offset + tail] = value;
            }
        }
        break;
    }
}

#undef LEO2_AVX2_TINY_GROUP_NOINLINE

void AVX2FF8LinearCombinationTiny(
    const void* const* sources,
    uint32_t source_count,
    void* const* destinations,
    const uint8_t* coefficient_logs,
    uint32_t output_count,
    uint64_t byte_count)
{
    if (!sources || source_count == 0 || source_count > 16 ||
        !destinations || !coefficient_logs || output_count == 0 ||
        output_count > 8)
        return;
    switch (output_count)
    {
    case 1:
        AVX2FF8EncodeOutputGroupTinyImpl<1>(
            sources, source_count, destinations,
            coefficient_logs, byte_count);
        break;
    case 2:
        AVX2FF8EncodeOutputGroupTinyImpl<2>(
            sources, source_count, destinations,
            coefficient_logs, byte_count);
        break;
    case 3:
        AVX2FF8EncodeOutputGroupTinyImpl<3>(
            sources, source_count, destinations,
            coefficient_logs, byte_count);
        break;
    case 4:
        AVX2FF8EncodeOutputGroupTinyImpl<4>(
            sources, source_count, destinations,
            coefficient_logs, byte_count);
        break;
    case 5:
        if (source_count == 5)
            AVX2FF8EncodeOutputGroupTinyImpl<5, 5>(
                sources, source_count, destinations,
                coefficient_logs, byte_count);
        else
            AVX2FF8EncodeOutputGroupTinyImpl<5>(
                sources, source_count, destinations,
                coefficient_logs, byte_count);
        break;
    case 6:
        AVX2FF8EncodeOutputGroupTinyImpl<6>(
            sources, source_count, destinations,
            coefficient_logs, byte_count);
        break;
    case 7:
        AVX2FF8EncodeOutputGroupTinyImpl<7>(
            sources, source_count, destinations,
            coefficient_logs, byte_count);
        break;
    case 8:
        AVX2FF8EncodeOutputGroupTinyImpl<8>(
            sources, source_count, destinations,
            coefficient_logs, byte_count);
        break;
    }
}

#endif

#if LEO2_EXPERIMENT_GENERAL_ONE_LOSS_DIRECT && \
    !defined(LEO2_AVX512_VARIANT) && !defined(LEO2_GFNI_VARIANT)
template<bool Add, bool Identity0, bool Identity1>
static void AVX2FF8LinearCombination2Loop(
    void* destination_pointer,
    const void* source0_pointer,
    const void* source1_pointer,
    uint16_t multiplier_log0,
    uint16_t multiplier_log1,
    uint64_t byte_count)
{
    LEO_DEBUG_ASSERT((byte_count & 31U) == 0);
    uint8_t* destination = static_cast<uint8_t*>(destination_pointer);
    const uint8_t* source0 = static_cast<const uint8_t*>(source0_pointer);
    const uint8_t* source1 = static_cast<const uint8_t*>(source1_pointer);
    const __m256i low0 = BroadcastTable(FF8Tables[multiplier_log0].low);
    const __m256i high0 = BroadcastTable(FF8Tables[multiplier_log0].high);
    const __m256i low1 = BroadcastTable(FF8Tables[multiplier_log1].low);
    const __m256i high1 = BroadcastTable(FF8Tables[multiplier_log1].high);
    for (uint64_t offset = 0; offset < byte_count; offset += 32)
    {
        const __m256i input0 = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(source0 + offset));
        const __m256i input1 = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(source1 + offset));
        const __m256i product0 = Identity0 ? input0 :
            AVX2FF8ProductVector(input0, low0, high0);
        const __m256i product1 = Identity1 ? input1 :
            AVX2FF8ProductVector(input1, low1, high1);
        __m256i result = _mm256_xor_si256(product0, product1);
        if (Add)
        {
            result = _mm256_xor_si256(result, _mm256_loadu_si256(
                reinterpret_cast<const __m256i*>(destination + offset)));
        }
        _mm256_storeu_si256(
            reinterpret_cast<__m256i*>(destination + offset), result);
    }
}

static LEO_FORCE_INLINE __m128i AVX2FF8LinearCombination2TailProduct(
    __m128i input,
    uint16_t multiplier_log,
    bool identity,
    __m128i nibble_mask)
{
    if (identity)
        return input;
    const FF8NibbleTable& table = FF8Tables[multiplier_log];
    const __m128i low_table = _mm_loadu_si128(
        reinterpret_cast<const __m128i*>(table.low));
    const __m128i high_table = _mm_loadu_si128(
        reinterpret_cast<const __m128i*>(table.high));
    const __m128i low_indices = _mm_and_si128(input, nibble_mask);
    const __m128i high_indices = _mm_and_si128(
        _mm_srli_epi64(input, 4), nibble_mask);
    return _mm_xor_si128(
        _mm_shuffle_epi8(low_table, low_indices),
        _mm_shuffle_epi8(high_table, high_indices));
}

template<bool Add>
static void AVX2FF8LinearCombination2Tail(
    uint8_t* destination,
    const uint8_t* source0,
    const uint8_t* source1,
    uint16_t multiplier_log0,
    uint16_t multiplier_log1,
    uint64_t byte_count)
{
    LEO_DEBUG_ASSERT(byte_count < 32);
    const bool identity0 = multiplier_log0 == 0 || multiplier_log0 == 255;
    const bool identity1 = multiplier_log1 == 0 || multiplier_log1 == 255;
    const __m128i nibble_mask = _mm_set1_epi8(0x0f);
    uint64_t offset = 0;
    if (byte_count >= 16)
    {
        const __m128i input0 = _mm_loadu_si128(
            reinterpret_cast<const __m128i*>(source0));
        const __m128i input1 = _mm_loadu_si128(
            reinterpret_cast<const __m128i*>(source1));
        __m128i result = _mm_xor_si128(
            AVX2FF8LinearCombination2TailProduct(
                input0, multiplier_log0, identity0, nibble_mask),
            AVX2FF8LinearCombination2TailProduct(
                input1, multiplier_log1, identity1, nibble_mask));
        if (Add)
        {
            result = _mm_xor_si128(result, _mm_loadu_si128(
                reinterpret_cast<const __m128i*>(destination)));
        }
        _mm_storeu_si128(reinterpret_cast<__m128i*>(destination), result);
        offset = 16;
    }
    if (byte_count - offset >= 8)
    {
        const __m128i input0 = _mm_loadl_epi64(
            reinterpret_cast<const __m128i*>(source0 + offset));
        const __m128i input1 = _mm_loadl_epi64(
            reinterpret_cast<const __m128i*>(source1 + offset));
        __m128i result = _mm_xor_si128(
            AVX2FF8LinearCombination2TailProduct(
                input0, multiplier_log0, identity0, nibble_mask),
            AVX2FF8LinearCombination2TailProduct(
                input1, multiplier_log1, identity1, nibble_mask));
        if (Add)
        {
            result = _mm_xor_si128(result, _mm_loadl_epi64(
                reinterpret_cast<const __m128i*>(destination + offset)));
        }
        _mm_storel_epi64(
            reinterpret_cast<__m128i*>(destination + offset), result);
        offset += 8;
    }
    while (offset < byte_count)
    {
        const uint8_t product0 = identity0 ? source0[offset] :
            FF8Product(multiplier_log0, source0[offset]);
        const uint8_t product1 = identity1 ? source1[offset] :
            FF8Product(multiplier_log1, source1[offset]);
        const uint8_t result = static_cast<uint8_t>(product0 ^ product1);
        if (Add)
            destination[offset] ^= result;
        else
            destination[offset] = result;
        ++offset;
    }
}

template<bool Add>
static void AVX2FF8LinearCombination2Mode(
    void* destination,
    const void* source0,
    const void* source1,
    uint16_t multiplier_log0,
    uint16_t multiplier_log1,
    uint64_t byte_count)
{
    static const uint16_t kSuppressed = UINT16_MAX;
    if (byte_count == 0)
        return;
    if (multiplier_log0 == kSuppressed && multiplier_log1 == kSuppressed)
    {
        if (!Add)
            memset(destination, 0, static_cast<size_t>(byte_count));
        return;
    }
    if (multiplier_log0 == kSuppressed || multiplier_log1 == kSuppressed)
    {
        const bool use_second = multiplier_log0 == kSuppressed;
        const void* source = use_second ? source1 : source0;
        const uint16_t multiplier_log = use_second
            ? multiplier_log1 : multiplier_log0;
        if (multiplier_log == 0 || multiplier_log == 255)
        {
            if (Add)
                AVX2XorMemory(destination, source, byte_count);
            else
                memcpy(destination, source, static_cast<size_t>(byte_count));
        }
        else if (Add)
        {
            AVX2FF8MultiplyAdd(
                destination, source, multiplier_log, byte_count);
        }
        else
        {
            AVX2FF8Multiply(destination, source, multiplier_log, byte_count);
        }
        return;
    }

    const bool identity0 = multiplier_log0 == 0 || multiplier_log0 == 255;
    const bool identity1 = multiplier_log1 == 0 || multiplier_log1 == 255;
    const uint64_t prefix_bytes = byte_count & ~UINT64_C(31);
    if (prefix_bytes != 0)
    {
        if (identity0)
        {
            if (identity1)
            {
                AVX2FF8LinearCombination2Loop<Add, true, true>(
                    destination, source0, source1,
                    multiplier_log0, multiplier_log1, prefix_bytes);
            }
            else
            {
                AVX2FF8LinearCombination2Loop<Add, true, false>(
                    destination, source0, source1,
                    multiplier_log0, multiplier_log1, prefix_bytes);
            }
        }
        else if (identity1)
        {
            AVX2FF8LinearCombination2Loop<Add, false, true>(
                destination, source0, source1,
                multiplier_log0, multiplier_log1, prefix_bytes);
        }
        else
        {
            AVX2FF8LinearCombination2Loop<Add, false, false>(
                destination, source0, source1,
                multiplier_log0, multiplier_log1, prefix_bytes);
        }
    }

    const uint64_t tail_bytes = byte_count - prefix_bytes;
    if (tail_bytes != 0)
    {
        AVX2FF8LinearCombination2Tail<Add>(
            static_cast<uint8_t*>(destination) + prefix_bytes,
            static_cast<const uint8_t*>(source0) + prefix_bytes,
            static_cast<const uint8_t*>(source1) + prefix_bytes,
            multiplier_log0, multiplier_log1, tail_bytes);
    }
}

static void AVX2FF8LinearCombination2(
    void* destination,
    const void* source0,
    const void* source1,
    uint16_t multiplier_log0,
    uint16_t multiplier_log1,
    bool add,
    uint64_t byte_count)
{
    if (add)
    {
        AVX2FF8LinearCombination2Mode<true>(
            destination, source0, source1,
            multiplier_log0, multiplier_log1, byte_count);
    }
    else
    {
        AVX2FF8LinearCombination2Mode<false>(
            destination, source0, source1,
            multiplier_log0, multiplier_log1, byte_count);
    }
}
#endif

#if LEO2_EXPERIMENT_GENERAL_ONE_LOSS_DIRECT && \
    !defined(LEO2_GFNI_VARIANT)
static LEO_FORCE_INLINE __m128i AVX2FF8FourTinyProduct128(
    __m128i input,
    __m128i low_table,
    __m128i high_table,
    __m128i nibble_mask)
{
    const __m128i low_indices = _mm_and_si128(input, nibble_mask);
    const __m128i high_indices = _mm_and_si128(
        _mm_srli_epi64(input, 4), nibble_mask);
    return _mm_xor_si128(
        _mm_shuffle_epi8(low_table, low_indices),
        _mm_shuffle_epi8(high_table, high_indices));
}

static LEO_FORCE_INLINE void AVX2FF8FourTinyConsume(__m128i& value)
{
#if (defined(__GNUC__) || defined(__clang__)) && !defined(_MSC_VER)
    // Bound temporary-product lifetimes while the nine reusable vector
    // constants are live.  This constraint emits no instruction.
    __asm__ __volatile__("" : "+x"(value));
#else
    (void)value;
#endif
}

static LEO_FORCE_INLINE uint8_t AVX2FF8FourTinyProduct(
    uint8_t input,
    uint16_t multiplier_log)
{
    // GF8 logarithms are modulo 255, so both representations name one.
    if (multiplier_log == 0 || multiplier_log == 255)
        return input;
    return FF8Product(multiplier_log, input);
}

template<bool Add>
static void AVX2FF8LinearCombination4TinyMode(
    void* destination_pointer,
    const void* const* source_pointers,
    const uint16_t* multiplier_logs,
    uint64_t byte_count)
{
    LEO_DEBUG_ASSERT(byte_count <= 63);
    uint8_t* const destination = static_cast<uint8_t*>(destination_pointer);
    const uint8_t* sources[4] = {
        static_cast<const uint8_t*>(source_pointers[0]),
        static_cast<const uint8_t*>(source_pointers[1]),
        static_cast<const uint8_t*>(source_pointers[2]),
        static_cast<const uint8_t*>(source_pointers[3])
    };
    uint64_t offset = 0;
    if (byte_count >= 8)
    {
        /*
            Keep all four table pairs and the shared mask live across the
            straight-line blocks.  The optional eight-byte tail is handled
            first so the compiler need not preserve those constants across a
            later tail branch.
        */
        const FF8NibbleTable& table0 = FF8Tables[multiplier_logs[0]];
        const FF8NibbleTable& table1 = FF8Tables[multiplier_logs[1]];
        const FF8NibbleTable& table2 = FF8Tables[multiplier_logs[2]];
        const FF8NibbleTable& table3 = FF8Tables[multiplier_logs[3]];
        const __m128i low0 = _mm_loadu_si128(
            reinterpret_cast<const __m128i*>(table0.low));
        const __m128i high0 = _mm_loadu_si128(
            reinterpret_cast<const __m128i*>(table0.high));
        const __m128i low1 = _mm_loadu_si128(
            reinterpret_cast<const __m128i*>(table1.low));
        const __m128i high1 = _mm_loadu_si128(
            reinterpret_cast<const __m128i*>(table1.high));
        const __m128i low2 = _mm_loadu_si128(
            reinterpret_cast<const __m128i*>(table2.low));
        const __m128i high2 = _mm_loadu_si128(
            reinterpret_cast<const __m128i*>(table2.high));
        const __m128i low3 = _mm_loadu_si128(
            reinterpret_cast<const __m128i*>(table3.low));
        const __m128i high3 = _mm_loadu_si128(
            reinterpret_cast<const __m128i*>(table3.high));
        __m128i nibble_mask = _mm_set1_epi8(0x0f);
#if (defined(__GNUC__) || defined(__clang__)) && !defined(_MSC_VER)
        __asm__ __volatile__("" : "+x"(nibble_mask));
#endif

        const uint64_t complete_bytes = byte_count & ~UINT64_C(15);
        const bool has_xmm_tail = byte_count - complete_bytes >= 8;
        if (has_xmm_tail)
        {
            const uint64_t tail_offset = complete_bytes;
            __m128i result = AVX2FF8FourTinyProduct128(
                _mm_loadl_epi64(reinterpret_cast<const __m128i*>(
                    sources[0] + tail_offset)),
                low0, high0, nibble_mask);
            AVX2FF8FourTinyConsume(result);
            result = _mm_xor_si128(result, AVX2FF8FourTinyProduct128(
                _mm_loadl_epi64(reinterpret_cast<const __m128i*>(
                    sources[1] + tail_offset)),
                low1, high1, nibble_mask));
            AVX2FF8FourTinyConsume(result);
            result = _mm_xor_si128(result, AVX2FF8FourTinyProduct128(
                _mm_loadl_epi64(reinterpret_cast<const __m128i*>(
                    sources[2] + tail_offset)),
                low2, high2, nibble_mask));
            AVX2FF8FourTinyConsume(result);
            result = _mm_xor_si128(result, AVX2FF8FourTinyProduct128(
                _mm_loadl_epi64(reinterpret_cast<const __m128i*>(
                    sources[3] + tail_offset)),
                low3, high3, nibble_mask));
            if (Add)
            {
                result = _mm_xor_si128(result, _mm_loadl_epi64(
                    reinterpret_cast<const __m128i*>(
                        destination + tail_offset)));
            }
            _mm_storel_epi64(reinterpret_cast<__m128i*>(
                destination + tail_offset), result);
        }

        while (offset < complete_bytes)
        {
            __m128i result = AVX2FF8FourTinyProduct128(
                _mm_loadu_si128(reinterpret_cast<const __m128i*>(
                    sources[0] + offset)),
                low0, high0, nibble_mask);
            AVX2FF8FourTinyConsume(result);
            result = _mm_xor_si128(result, AVX2FF8FourTinyProduct128(
                _mm_loadu_si128(reinterpret_cast<const __m128i*>(
                    sources[1] + offset)),
                low1, high1, nibble_mask));
            AVX2FF8FourTinyConsume(result);
            result = _mm_xor_si128(result, AVX2FF8FourTinyProduct128(
                _mm_loadu_si128(reinterpret_cast<const __m128i*>(
                    sources[2] + offset)),
                low2, high2, nibble_mask));
            AVX2FF8FourTinyConsume(result);
            result = _mm_xor_si128(result, AVX2FF8FourTinyProduct128(
                _mm_loadu_si128(reinterpret_cast<const __m128i*>(
                    sources[3] + offset)),
                low3, high3, nibble_mask));
            if (Add)
            {
                result = _mm_xor_si128(result, _mm_loadu_si128(
                    reinterpret_cast<const __m128i*>(
                        destination + offset)));
            }
            _mm_storeu_si128(reinterpret_cast<__m128i*>(
                destination + offset), result);
            offset += 16;
        }
        if (has_xmm_tail)
            offset += 8;
    }

    while (offset < byte_count)
    {
        uint8_t result = AVX2FF8FourTinyProduct(
            sources[0][offset], multiplier_logs[0]);
        for (unsigned source = 1; source < 4; ++source)
        {
            result ^= AVX2FF8FourTinyProduct(
                sources[source][offset], multiplier_logs[source]);
        }
        if (Add)
            destination[offset] ^= result;
        else
            destination[offset] = result;
        ++offset;
    }
}

static void AVX2FF8LinearCombination4Tiny(
    void* destination,
    const void* const* sources,
    const uint16_t* multiplier_logs,
    bool add,
    uint64_t byte_count)
{
    if (byte_count == 0)
        return;
    LEO_DEBUG_ASSERT(destination && sources && multiplier_logs);
    LEO_DEBUG_ASSERT(byte_count <= 63);
    if (add)
    {
        AVX2FF8LinearCombination4TinyMode<true>(
            destination, sources, multiplier_logs, byte_count);
    }
    else
    {
        AVX2FF8LinearCombination4TinyMode<false>(
            destination, sources, multiplier_logs, byte_count);
    }
}
#endif // LEO2_EXPERIMENT_GENERAL_ONE_LOSS_DIRECT && !LEO2_GFNI_VARIANT

static void AVX2FF8MultiplyAdd2Sources2Outputs(
    void* const* destination_pointers,
    const void* source0_pointer,
    const void* source1_pointer,
    const uint16_t* multiplier_logs0,
    const uint16_t* multiplier_logs1,
    uint64_t byte_count)
{
    uint8_t* destination0 = static_cast<uint8_t*>(destination_pointers[0]);
    uint8_t* destination1 = static_cast<uint8_t*>(destination_pointers[1]);
    const uint8_t* source0 = static_cast<const uint8_t*>(source0_pointer);
    const uint8_t* source1 = static_cast<const uint8_t*>(source1_pointer);
    const __m256i low00 = BroadcastTable(FF8Tables[multiplier_logs0[0]].low);
    const __m256i high00 = BroadcastTable(FF8Tables[multiplier_logs0[0]].high);
    const __m256i low01 = BroadcastTable(FF8Tables[multiplier_logs0[1]].low);
    const __m256i high01 = BroadcastTable(FF8Tables[multiplier_logs0[1]].high);
    const __m256i low10 = BroadcastTable(FF8Tables[multiplier_logs1[0]].low);
    const __m256i high10 = BroadcastTable(FF8Tables[multiplier_logs1[0]].high);
    const __m256i low11 = BroadcastTable(FF8Tables[multiplier_logs1[1]].low);
    const __m256i high11 = BroadcastTable(FF8Tables[multiplier_logs1[1]].high);
    uint64_t offset = 0;
    while (byte_count - offset >= 32)
    {
        const __m256i input0 = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(source0 + offset));
        const __m256i input1 = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(source1 + offset));
        __m256i result = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(destination0 + offset));
        result = _mm256_xor_si256(result,
            AVX2FF8ProductVector(input0, low00, high00));
        result = _mm256_xor_si256(result,
            AVX2FF8ProductVector(input1, low10, high10));
        _mm256_storeu_si256(
            reinterpret_cast<__m256i*>(destination0 + offset), result);
        result = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(destination1 + offset));
        result = _mm256_xor_si256(result,
            AVX2FF8ProductVector(input0, low01, high01));
        result = _mm256_xor_si256(result,
            AVX2FF8ProductVector(input1, low11, high11));
        _mm256_storeu_si256(
            reinterpret_cast<__m256i*>(destination1 + offset), result);
        offset += 32;
    }
    while (offset < byte_count)
    {
        destination0[offset] ^=
            FF8Product(multiplier_logs0[0], source0[offset]) ^
            FF8Product(multiplier_logs1[0], source1[offset]);
        destination1[offset] ^=
            FF8Product(multiplier_logs0[1], source0[offset]) ^
            FF8Product(multiplier_logs1[1], source1[offset]);
        ++offset;
    }
}

#endif


#if defined(LEO2_GFNI_VARIANT) && defined(LEO_HAS_FF8)
/*
    One out-of-place inverse radix-eight group: three transform layers per
    load/store round instead of two.

    The staging kernel this replaces owns about 40% of a legacy-high GF8 encode
    call and runs within a few percent of a kernel that only moves the bytes, so
    it is bound by load/store throughput, not arithmetic.  The way past a
    movement floor is to move less: eight loads and eight stores carry
    twenty-four shard-layers here against sixteen for two radix-four rounds, a
    third less traffic per unit of transform work.

    Register pressure is not the obstacle it looks like.  Out-of-place means the
    eight inputs are consumed into registers before any output is written, so
    eight data vectors are live rather than sixteen, and every multiplier matrix
    is a memory operand, which costs nothing on a memory-bound kernel.  This is
    why the same idea does not pay for the in-place ALU-bound kernels.
*/
static LEO_FORCE_INLINE void AVX2FF8Radix8Butterfly(
    __m256i& x, __m256i& y, uint16_t log)
{
    static const uint16_t kZeroSkew = 255;
    y = _mm256_xor_si256(y, x);
    if (log != kZeroSkew)
        x = _mm256_xor_si256(x, _mm256_gf2p8affine_epi64_epi8(
            y, _mm256_loadu_si256(reinterpret_cast<const __m256i*>(
                FF8Tables[log].low)), 0));
}

static LEO_FORCE_INLINE void AVX2FF8Radix8ForwardButterfly(
    __m256i& x, __m256i& y, uint16_t log)
{
    static const uint16_t kZeroSkew = 255;
    if (log != kZeroSkew)
        x = _mm256_xor_si256(x, _mm256_gf2p8affine_epi64_epi8(
            y, _mm256_loadu_si256(reinterpret_cast<const __m256i*>(
                FF8Tables[log].low)), 0));
    y = _mm256_xor_si256(y, x);
}

#define LEO2_R8_LOAD(NAME, LANE) \
    __m256i NAME = _mm256_loadu_si256( \
        reinterpret_cast<const __m256i*>(in[LANE] + offset));
#define LEO2_R8_STORE(LANE, NAME) \
    _mm256_storeu_si256(reinterpret_cast<__m256i*>(out[LANE] + offset), NAME);

static void AVX2FF8IFFTButterfly8Out(
    const void* const* inputs,
    void* const* outputs,
    const uint8_t* skews,
    uint64_t byte_count)
{
    const uint8_t* in[8];
    uint8_t* out[8];
    for (unsigned lane = 0; lane < 8; ++lane)
    {
        in[lane] = static_cast<const uint8_t*>(inputs[lane]);
        out[lane] = static_cast<uint8_t*>(outputs[lane]);
    }
    const uint16_t s01 = skews[0], s23 = skews[1], s45 = skews[2];
    const uint16_t s67 = skews[3], s02 = skews[4], s46 = skews[5];
    const uint16_t s04 = skews[6];
    uint64_t offset = 0;
    while (byte_count - offset >= 32)
    {
        LEO2_R8_LOAD(x0, 0) LEO2_R8_LOAD(x1, 1)
        LEO2_R8_LOAD(x2, 2) LEO2_R8_LOAD(x3, 3)
        LEO2_R8_LOAD(x4, 4) LEO2_R8_LOAD(x5, 5)
        LEO2_R8_LOAD(x6, 6) LEO2_R8_LOAD(x7, 7)
        AVX2FF8Radix8Butterfly(x0, x1, s01);
        AVX2FF8Radix8Butterfly(x2, x3, s23);
        AVX2FF8Radix8Butterfly(x4, x5, s45);
        AVX2FF8Radix8Butterfly(x6, x7, s67);
        AVX2FF8Radix8Butterfly(x0, x2, s02);
        AVX2FF8Radix8Butterfly(x1, x3, s02);
        AVX2FF8Radix8Butterfly(x4, x6, s46);
        AVX2FF8Radix8Butterfly(x5, x7, s46);
        AVX2FF8Radix8Butterfly(x0, x4, s04);
        AVX2FF8Radix8Butterfly(x1, x5, s04);
        AVX2FF8Radix8Butterfly(x2, x6, s04);
        AVX2FF8Radix8Butterfly(x3, x7, s04);
        LEO2_R8_STORE(0, x0) LEO2_R8_STORE(1, x1)
        LEO2_R8_STORE(2, x2) LEO2_R8_STORE(3, x3)
        LEO2_R8_STORE(4, x4) LEO2_R8_STORE(5, x5)
        LEO2_R8_STORE(6, x6) LEO2_R8_STORE(7, x7)
        offset += 32;
    }
    while (offset < byte_count)
    {
        uint8_t v[8];
        for (unsigned lane = 0; lane < 8; ++lane)
            v[lane] = in[lane][offset];
        static const uint16_t kZeroSkewScalar = 255;
        static const unsigned pairs[12][3] = {
            {0,1,0},{2,3,1},{4,5,2},{6,7,3},
            {0,2,4},{1,3,4},{4,6,5},{5,7,5},
            {0,4,6},{1,5,6},{2,6,6},{3,7,6}
        };
        for (unsigned p = 0; p < 12; ++p)
        {
            const unsigned a = pairs[p][0], b = pairs[p][1];
            const uint16_t log = skews[pairs[p][2]];
            v[b] ^= v[a];
            if (log != kZeroSkewScalar)
                v[a] ^= FF8Product(log, v[b]);
        }
        for (unsigned lane = 0; lane < 8; ++lane)
            out[lane][offset] = v[lane];
        ++offset;
    }
}

/*
    Forward mirror.  The group order reverses — largest distance first — and the
    butterfly becomes x ^= mul(y, log); y ^= x, both taken from the forward
    branch of AVX2FF8Butterfly4Out, which applies the distance-2d skew before
    the two distance-d skews.  Skews arrive largest distance first:
    skew[4d], skew[2d], skew[6d], skew[d], skew[3d], skew[5d], skew[7d].
*/
static void AVX2FF8FFTButterfly8Out(
    const void* const* inputs,
    void* const* outputs,
    const uint8_t* skews,
    uint64_t byte_count)
{
    const uint8_t* in[8];
    uint8_t* out[8];
    for (unsigned lane = 0; lane < 8; ++lane)
    {
        in[lane] = static_cast<const uint8_t*>(inputs[lane]);
        out[lane] = static_cast<uint8_t*>(outputs[lane]);
    }
    const uint16_t s04 = skews[0];
    const uint16_t s02 = skews[1], s46 = skews[2];
    const uint16_t s01 = skews[3], s23 = skews[4];
    const uint16_t s45 = skews[5], s67 = skews[6];
    uint64_t offset = 0;
    while (byte_count - offset >= 32)
    {
        LEO2_R8_LOAD(x0, 0) LEO2_R8_LOAD(x1, 1)
        LEO2_R8_LOAD(x2, 2) LEO2_R8_LOAD(x3, 3)
        LEO2_R8_LOAD(x4, 4) LEO2_R8_LOAD(x5, 5)
        LEO2_R8_LOAD(x6, 6) LEO2_R8_LOAD(x7, 7)
        AVX2FF8Radix8ForwardButterfly(x0, x4, s04);
        AVX2FF8Radix8ForwardButterfly(x1, x5, s04);
        AVX2FF8Radix8ForwardButterfly(x2, x6, s04);
        AVX2FF8Radix8ForwardButterfly(x3, x7, s04);
        AVX2FF8Radix8ForwardButterfly(x0, x2, s02);
        AVX2FF8Radix8ForwardButterfly(x1, x3, s02);
        AVX2FF8Radix8ForwardButterfly(x4, x6, s46);
        AVX2FF8Radix8ForwardButterfly(x5, x7, s46);
        AVX2FF8Radix8ForwardButterfly(x0, x1, s01);
        AVX2FF8Radix8ForwardButterfly(x2, x3, s23);
        AVX2FF8Radix8ForwardButterfly(x4, x5, s45);
        AVX2FF8Radix8ForwardButterfly(x6, x7, s67);
        LEO2_R8_STORE(0, x0) LEO2_R8_STORE(1, x1)
        LEO2_R8_STORE(2, x2) LEO2_R8_STORE(3, x3)
        LEO2_R8_STORE(4, x4) LEO2_R8_STORE(5, x5)
        LEO2_R8_STORE(6, x6) LEO2_R8_STORE(7, x7)
        offset += 32;
    }
    while (offset < byte_count)
    {
        uint8_t v[8];
        for (unsigned lane = 0; lane < 8; ++lane)
            v[lane] = in[lane][offset];
        static const uint16_t kZeroSkewScalar = 255;
        static const unsigned pairs[12][3] = {
            {0,4,0},{1,5,0},{2,6,0},{3,7,0},
            {0,2,1},{1,3,1},{4,6,2},{5,7,2},
            {0,1,3},{2,3,4},{4,5,5},{6,7,6}
        };
        for (unsigned p = 0; p < 12; ++p)
        {
            const unsigned a = pairs[p][0], b = pairs[p][1];
            const uint16_t log = skews[pairs[p][2]];
            if (log != kZeroSkewScalar)
                v[a] ^= FF8Product(log, v[b]);
            v[b] ^= v[a];
        }
        for (unsigned lane = 0; lane < 8; ++lane)
            out[lane][offset] = v[lane];
        ++offset;
    }
}

#undef LEO2_R8_LOAD
#undef LEO2_R8_STORE
#endif

#if defined(LEO_HAS_FF8) && \
    !defined(LEO2_DIAGNOSTIC_DISABLE_HIGH_T8_VECTOR) && \
    !defined(LEO2_AVX512_VARIANT) && !defined(LEO2_GFNI_VARIANT)

/*
    The AVX-512 implementation below keeps two 32-byte vectors for every T=8
    coordinate live at once.  That is appropriate with 32 architectural vector
    registers but spills on ordinary AVX2.  The regular AVX2 callback processes
    one 32-byte slice at a time so the eight code symbols, two nibble tables,
    and product temporaries fit in the 16-register AVX2 file.
*/
static LEO_FORCE_INLINE void AVX2FF8T8VectorXor(
    __m256i& destination,
    __m256i source)
{
    destination = _mm256_xor_si256(destination, source);
}

static LEO_FORCE_INLINE void AVX2FF8T8VectorMulAddPrepared(
    __m256i& destination,
    __m256i source,
    __m256i low_table,
    __m256i high_table)
{
    destination = _mm256_xor_si256(destination,
        AVX2FF8ProductVector(source, low_table, high_table));
}

static LEO_FORCE_INLINE void AVX2FF8T8VectorMulAdd(
    __m256i& destination,
    __m256i source,
    uint16_t log)
{
    const FF8NibbleTable& table = FF8Tables[log];
    AVX2FF8T8VectorMulAddPrepared(destination, source,
        BroadcastTable(table.low), BroadcastTable(table.high));
}

static LEO_FORCE_INLINE void AVX2FF8T8VectorIFFTRadix4(
    __m256i& value0,
    __m256i& value1,
    __m256i& value2,
    __m256i& value3,
    uint16_t log01,
    uint16_t log23,
    uint16_t log02)
{
    static const uint16_t kZeroSkew = 255;
    AVX2FF8T8VectorXor(value1, value0);
    if (log01 != kZeroSkew)
        AVX2FF8T8VectorMulAdd(value0, value1, log01);
    AVX2FF8T8VectorXor(value3, value2);
    if (log23 != kZeroSkew)
        AVX2FF8T8VectorMulAdd(value2, value3, log23);
    AVX2FF8T8VectorXor(value2, value0);
    AVX2FF8T8VectorXor(value3, value1);
    if (log02 != kZeroSkew)
    {
        const FF8NibbleTable& table = FF8Tables[log02];
        const __m256i low_table = BroadcastTable(table.low);
        const __m256i high_table = BroadcastTable(table.high);
        AVX2FF8T8VectorMulAddPrepared(
            value0, value2, low_table, high_table);
        AVX2FF8T8VectorMulAddPrepared(
            value1, value3, low_table, high_table);
    }
}

static LEO_FORCE_INLINE void AVX2FF8T8VectorIFFTDistance4(
    __m256i& value0,
    __m256i& value1,
    __m256i& value2,
    __m256i& value3,
    __m256i& value4,
    __m256i& value5,
    __m256i& value6,
    __m256i& value7,
    uint16_t log)
{
    static const uint16_t kZeroSkew = 255;
    AVX2FF8T8VectorXor(value4, value0);
    AVX2FF8T8VectorXor(value5, value1);
    AVX2FF8T8VectorXor(value6, value2);
    AVX2FF8T8VectorXor(value7, value3);
    if (log == kZeroSkew)
        return;
    const FF8NibbleTable& table = FF8Tables[log];
    const __m256i low_table = BroadcastTable(table.low);
    const __m256i high_table = BroadcastTable(table.high);
    AVX2FF8T8VectorMulAddPrepared(value0, value4, low_table, high_table);
    AVX2FF8T8VectorMulAddPrepared(value1, value5, low_table, high_table);
    AVX2FF8T8VectorMulAddPrepared(value2, value6, low_table, high_table);
    AVX2FF8T8VectorMulAddPrepared(value3, value7, low_table, high_table);
}

static LEO_FORCE_INLINE void AVX2FF8T8VectorFFTRadix4Distance2(
    __m256i& value0,
    __m256i& value1,
    __m256i& value2,
    __m256i& value3,
    __m256i& value4,
    __m256i& value5,
    __m256i& value6,
    __m256i& value7,
    uint16_t log01,
    uint16_t log23,
    uint16_t log02)
{
    static const uint16_t kZeroSkew = 255;
    if (log02 != kZeroSkew)
    {
        const FF8NibbleTable& table = FF8Tables[log02];
        const __m256i low_table = BroadcastTable(table.low);
        const __m256i high_table = BroadcastTable(table.high);
        AVX2FF8T8VectorMulAddPrepared(
            value0, value4, low_table, high_table);
        AVX2FF8T8VectorMulAddPrepared(
            value1, value5, low_table, high_table);
        AVX2FF8T8VectorMulAddPrepared(
            value2, value6, low_table, high_table);
        AVX2FF8T8VectorMulAddPrepared(
            value3, value7, low_table, high_table);
    }
    AVX2FF8T8VectorXor(value4, value0);
    AVX2FF8T8VectorXor(value5, value1);
    AVX2FF8T8VectorXor(value6, value2);
    AVX2FF8T8VectorXor(value7, value3);
    if (log01 != kZeroSkew)
    {
        const FF8NibbleTable& table = FF8Tables[log01];
        const __m256i low_table = BroadcastTable(table.low);
        const __m256i high_table = BroadcastTable(table.high);
        AVX2FF8T8VectorMulAddPrepared(
            value0, value2, low_table, high_table);
        AVX2FF8T8VectorMulAddPrepared(
            value1, value3, low_table, high_table);
    }
    AVX2FF8T8VectorXor(value2, value0);
    AVX2FF8T8VectorXor(value3, value1);
    if (log23 != kZeroSkew)
    {
        const FF8NibbleTable& table = FF8Tables[log23];
        const __m256i low_table = BroadcastTable(table.low);
        const __m256i high_table = BroadcastTable(table.high);
        AVX2FF8T8VectorMulAddPrepared(
            value4, value6, low_table, high_table);
        AVX2FF8T8VectorMulAddPrepared(
            value5, value7, low_table, high_table);
    }
    AVX2FF8T8VectorXor(value6, value4);
    AVX2FF8T8VectorXor(value7, value5);
}

static LEO_FORCE_INLINE void AVX2FF8T8VectorFFT2(
    __m256i& value0,
    __m256i& value1,
    uint16_t log)
{
    if (log != 255)
        AVX2FF8T8VectorMulAdd(value0, value1, log);
    AVX2FF8T8VectorXor(value1, value0);
}

static LEO_FORCE_INLINE void AVX2FF8T8VectorScalarIFFT2(
    uint8_t& value0,
    uint8_t& value1,
    uint16_t log)
{
    value1 ^= value0;
    if (log != 255)
        value0 ^= FF8Product(log, value1);
}

static LEO_FORCE_INLINE void AVX2FF8T8VectorScalarFFT2(
    uint8_t& value0,
    uint8_t& value1,
    uint16_t log)
{
    if (log != 255)
        value0 ^= FF8Product(log, value1);
    value1 ^= value0;
}

static void AVX2FF8HighEncodeT8Vector(
    const void* const* data,
    void* const* work,
    bool shortened,
    const uint8_t* inverse_skew,
    const uint8_t* forward_skew,
    uint64_t byte_count)
{
    alignas(32) static const uint8_t kZeroShard[32] = {};
    const uint8_t* input[8];
    uint8_t* output[8];
    for (unsigned lane = 0; lane < 7; ++lane)
    {
        input[lane] = static_cast<const uint8_t*>(data[lane]);
        output[lane] = static_cast<uint8_t*>(work[lane]);
    }
    input[7] = shortened
        ? kZeroShard : static_cast<const uint8_t*>(data[7]);
    output[7] = static_cast<uint8_t*>(work[7]);
    const uint64_t lane7_offset_mask = shortened ? 0 : ~uint64_t(0);

    uint64_t offset = 0;
    while (byte_count - offset >= 32)
    {
        __m256i value0 = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(input[0] + offset));
        __m256i value1 = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(input[1] + offset));
        __m256i value2 = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(input[2] + offset));
        __m256i value3 = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(input[3] + offset));
        __m256i value4 = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(input[4] + offset));
        __m256i value5 = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(input[5] + offset));
        __m256i value6 = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(input[6] + offset));
        const uint64_t lane7_offset = offset & lane7_offset_mask;
        __m256i value7 = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(
                input[7] + lane7_offset));

        AVX2FF8T8VectorIFFTRadix4(
            value0, value1, value2, value3,
            inverse_skew[1], inverse_skew[3], inverse_skew[2]);
        AVX2FF8T8VectorIFFTRadix4(
            value4, value5, value6, value7,
            inverse_skew[5], inverse_skew[7], inverse_skew[6]);
        AVX2FF8T8VectorIFFTDistance4(
            value0, value1, value2, value3,
            value4, value5, value6, value7, inverse_skew[4]);

        AVX2FF8T8VectorFFTRadix4Distance2(
            value0, value1, value2, value3,
            value4, value5, value6, value7,
            forward_skew[2], forward_skew[6], forward_skew[4]);
        AVX2FF8T8VectorFFT2(value0, value1, forward_skew[1]);
        AVX2FF8T8VectorFFT2(value2, value3, forward_skew[3]);
        AVX2FF8T8VectorFFT2(value4, value5, forward_skew[5]);
        AVX2FF8T8VectorFFT2(value6, value7, forward_skew[7]);

        _mm256_storeu_si256(
            reinterpret_cast<__m256i*>(output[0] + offset), value0);
        _mm256_storeu_si256(
            reinterpret_cast<__m256i*>(output[1] + offset), value1);
        _mm256_storeu_si256(
            reinterpret_cast<__m256i*>(output[2] + offset), value2);
        _mm256_storeu_si256(
            reinterpret_cast<__m256i*>(output[3] + offset), value3);
        _mm256_storeu_si256(
            reinterpret_cast<__m256i*>(output[4] + offset), value4);
        _mm256_storeu_si256(
            reinterpret_cast<__m256i*>(output[5] + offset), value5);
        _mm256_storeu_si256(
            reinterpret_cast<__m256i*>(output[6] + offset), value6);
        _mm256_storeu_si256(
            reinterpret_cast<__m256i*>(output[7] + offset), value7);
        offset += 32;
    }

    const size_t tail_bytes = static_cast<size_t>(byte_count - offset);
    if (tail_bytes != 0)
    {
        /*
            When at least one complete vector exists, evaluate the final
            in-range 32-byte window.  Its overlap with the preceding vector is
            harmless because every byte position is independent, and this
            avoids all tail copies for byte counts above 32.  A genuinely
            short shard still uses one zero-padded staging vector.
        */
        alignas(32) uint8_t tail_input[8][32] = {};
        alignas(32) uint8_t tail_output[8][32];
        const uint8_t* tail_source[8];
        uint8_t* tail_destination[8];
        const bool staged_tail = offset == 0;
        if (staged_tail)
        {
            for (unsigned lane = 0; lane < 7; ++lane)
            {
                std::memcpy(
                    tail_input[lane], input[lane], tail_bytes);
                tail_source[lane] = tail_input[lane];
                tail_destination[lane] = tail_output[lane];
            }
            if (!shortened)
                std::memcpy(
                    tail_input[7], input[7], tail_bytes);
            tail_source[7] = tail_input[7];
            tail_destination[7] = tail_output[7];
        }
        else
        {
            const uint64_t final_offset = byte_count - 32;
            for (unsigned lane = 0; lane < 7; ++lane)
            {
                tail_source[lane] = input[lane] + final_offset;
                tail_destination[lane] = output[lane] + final_offset;
            }
            tail_source[7] = input[7] +
                (shortened ? 0 : final_offset);
            tail_destination[7] = output[7] + final_offset;
        }

        __m256i value0 = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(tail_source[0]));
        __m256i value1 = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(tail_source[1]));
        __m256i value2 = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(tail_source[2]));
        __m256i value3 = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(tail_source[3]));
        __m256i value4 = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(tail_source[4]));
        __m256i value5 = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(tail_source[5]));
        __m256i value6 = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(tail_source[6]));
        __m256i value7 = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(tail_source[7]));

        AVX2FF8T8VectorIFFTRadix4(
            value0, value1, value2, value3,
            inverse_skew[1], inverse_skew[3], inverse_skew[2]);
        AVX2FF8T8VectorIFFTRadix4(
            value4, value5, value6, value7,
            inverse_skew[5], inverse_skew[7], inverse_skew[6]);
        AVX2FF8T8VectorIFFTDistance4(
            value0, value1, value2, value3,
            value4, value5, value6, value7, inverse_skew[4]);

        AVX2FF8T8VectorFFTRadix4Distance2(
            value0, value1, value2, value3,
            value4, value5, value6, value7,
            forward_skew[2], forward_skew[6], forward_skew[4]);
        AVX2FF8T8VectorFFT2(value0, value1, forward_skew[1]);
        AVX2FF8T8VectorFFT2(value2, value3, forward_skew[3]);
        AVX2FF8T8VectorFFT2(value4, value5, forward_skew[5]);
        AVX2FF8T8VectorFFT2(value6, value7, forward_skew[7]);

        _mm256_storeu_si256(
            reinterpret_cast<__m256i*>(tail_destination[0]), value0);
        _mm256_storeu_si256(
            reinterpret_cast<__m256i*>(tail_destination[1]), value1);
        _mm256_storeu_si256(
            reinterpret_cast<__m256i*>(tail_destination[2]), value2);
        _mm256_storeu_si256(
            reinterpret_cast<__m256i*>(tail_destination[3]), value3);
        _mm256_storeu_si256(
            reinterpret_cast<__m256i*>(tail_destination[4]), value4);
        _mm256_storeu_si256(
            reinterpret_cast<__m256i*>(tail_destination[5]), value5);
        _mm256_storeu_si256(
            reinterpret_cast<__m256i*>(tail_destination[6]), value6);
        _mm256_storeu_si256(
            reinterpret_cast<__m256i*>(tail_destination[7]), value7);
        if (staged_tail)
            for (unsigned lane = 0; lane < 8; ++lane)
                std::memcpy(
                    output[lane], tail_output[lane], tail_bytes);
    }
}

#if defined(_MSC_VER)
# define LEO2_AVX2_T8_ENTRY __declspec(noinline)
#elif defined(__GNUC__) || defined(__clang__)
# define LEO2_AVX2_T8_ENTRY \
    __attribute__((noinline, section(".text.leo2_t8")))
#else
# define LEO2_AVX2_T8_ENTRY
#endif

#if LEO2_EXPERIMENT_HIGH_T8_TWO_BLOCK_BINDING
/*
    Preserve the already-qualified exact K=9/R=5 body and section placement.
    R=6..8 below use the same first-block transform with additional generator
    rows, but must not perturb this promoted terminal's register allocation.
*/
static LEO2_AVX2_T8_ENTRY void AVX2FF8HighEncodeK9R5T8Multiple32(
    const void* const* data,
    void* const* work,
    uint64_t byte_count)
{
    static const uint16_t kInverse1 = 153;
    static const uint16_t kInverse2 = 17;
    static const uint16_t kInverse3 = 102;
    static const uint16_t kInverse4 = 85;
    static const uint16_t kInverse5 = 51;
    static const uint16_t kInverse6 = 34;
    static const uint16_t kInverse7 = 187;
    static const uint16_t kForward1 = 255;
    static const uint16_t kForward2 = 255;
    static const uint16_t kForward3 = 85;
    static const uint16_t kForward4 = 255;
    static const uint16_t kForward5 = 17;
    static const uint16_t kForward6 = 85;
    static const uint16_t kForward7 = 34;
    static const uint16_t kTail0 = 121;
    static const uint16_t kTail1 = 151;
    static const uint16_t kTail2 = 78;
    static const uint16_t kTail3 = 228;
    static const uint16_t kTail4 = 229;
    LEO_DEBUG_ASSERT(byte_count != 0 && (byte_count & 31U) == 0);
    if (byte_count == 0 || (byte_count & 31U) != 0)
        return;

    const uint8_t* input[9];
    uint8_t* output[5];
    for (unsigned lane = 0; lane < 9; ++lane)
        input[lane] = static_cast<const uint8_t*>(data[lane]);
    for (unsigned lane = 0; lane < 5; ++lane)
        output[lane] = static_cast<uint8_t*>(work[lane]);

    for (uint64_t offset = 0; offset < byte_count; offset += 32)
    {
        __m256i value0 = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(input[0] + offset));
        __m256i value1 = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(input[1] + offset));
        __m256i value2 = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(input[2] + offset));
        __m256i value3 = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(input[3] + offset));
        __m256i value4 = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(input[4] + offset));
        __m256i value5 = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(input[5] + offset));
        __m256i value6 = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(input[6] + offset));
        __m256i value7 = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(input[7] + offset));
        const __m256i tail = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(input[8] + offset));

        AVX2FF8T8VectorIFFTRadix4(
            value0, value1, value2, value3,
            kInverse1, kInverse3, kInverse2);
        AVX2FF8T8VectorIFFTRadix4(
            value4, value5, value6, value7,
            kInverse5, kInverse7, kInverse6);
        AVX2FF8T8VectorIFFTDistance4(
            value0, value1, value2, value3,
            value4, value5, value6, value7, kInverse4);

        AVX2FF8T8VectorFFTRadix4Distance2(
            value0, value1, value2, value3,
            value4, value5, value6, value7,
            kForward2, kForward6, kForward4);
        AVX2FF8T8VectorFFT2(value0, value1, kForward1);
        AVX2FF8T8VectorFFT2(value2, value3, kForward3);
        AVX2FF8T8VectorFFT2(value4, value5, kForward5);
        AVX2FF8T8VectorFFT2(value6, value7, kForward7);

        AVX2FF8T8VectorMulAdd(value0, tail, kTail0);
        AVX2FF8T8VectorMulAdd(value1, tail, kTail1);
        AVX2FF8T8VectorMulAdd(value2, tail, kTail2);
        AVX2FF8T8VectorMulAdd(value3, tail, kTail3);
        AVX2FF8T8VectorMulAdd(value4, tail, kTail4);

        _mm256_storeu_si256(
            reinterpret_cast<__m256i*>(output[0] + offset), value0);
        _mm256_storeu_si256(
            reinterpret_cast<__m256i*>(output[1] + offset), value1);
        _mm256_storeu_si256(
            reinterpret_cast<__m256i*>(output[2] + offset), value2);
        _mm256_storeu_si256(
            reinterpret_cast<__m256i*>(output[3] + offset), value3);
        _mm256_storeu_si256(
            reinterpret_cast<__m256i*>(output[4] + offset), value4);
    }
}

/*
    Exact K=9/R=6..8 legacy-high tiles.  The first eight systematic shards are a
    complete T=8 inverse-transform block.  The ninth shard is coordinate 16
    in the [32,24] parent and contributes the fixed systematic-generator
    column whose logarithms are listed below.  Applying that column after the
    first block's parity transform is algebraically identical to the second
    shifted inverse transform, but removes seven zero-row loads, the complete
    second IFFT, eight intermediate stores/loads, and unused punctured stores.
*/
template<unsigned OutputCount>
static LEO2_AVX2_T8_ENTRY void AVX2FF8HighEncodeK9T8Multiple32(
    const void* const* data,
    void* const* work,
    uint64_t byte_count)
{
    static const uint16_t kInverse1 = 153;
    static const uint16_t kInverse2 = 17;
    static const uint16_t kInverse3 = 102;
    static const uint16_t kInverse4 = 85;
    static const uint16_t kInverse5 = 51;
    static const uint16_t kInverse6 = 34;
    static const uint16_t kInverse7 = 187;
    static const uint16_t kForward1 = 255;
    static const uint16_t kForward2 = 255;
    static const uint16_t kForward3 = 85;
    static const uint16_t kForward4 = 255;
    static const uint16_t kForward5 = 17;
    static const uint16_t kForward6 = 85;
    static const uint16_t kForward7 = 34;
    static const uint16_t kTail0 = 121;
    static const uint16_t kTail1 = 151;
    static const uint16_t kTail2 = 78;
    static const uint16_t kTail3 = 228;
    static const uint16_t kTail4 = 229;
    static const uint16_t kTail5 = 94;
    static const uint16_t kTail6 = 57;
    static const uint16_t kTail7 = 147;
    static_assert(OutputCount >= 6 && OutputCount <= 8,
        "K9 extended terminal output count must be between six and eight");
    LEO_DEBUG_ASSERT(byte_count != 0 && (byte_count & 31U) == 0);
    if (byte_count == 0 || (byte_count & 31U) != 0)
        return;

    const uint8_t* input[9];
    uint8_t* output[OutputCount];
    for (unsigned lane = 0; lane < 9; ++lane)
        input[lane] = static_cast<const uint8_t*>(data[lane]);
    for (unsigned lane = 0; lane < OutputCount; ++lane)
        output[lane] = static_cast<uint8_t*>(work[lane]);

    for (uint64_t offset = 0; offset < byte_count; offset += 32)
    {
        __m256i value0 = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(input[0] + offset));
        __m256i value1 = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(input[1] + offset));
        __m256i value2 = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(input[2] + offset));
        __m256i value3 = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(input[3] + offset));
        __m256i value4 = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(input[4] + offset));
        __m256i value5 = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(input[5] + offset));
        __m256i value6 = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(input[6] + offset));
        __m256i value7 = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(input[7] + offset));
        const __m256i tail = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(input[8] + offset));

        AVX2FF8T8VectorIFFTRadix4(
            value0, value1, value2, value3,
            kInverse1, kInverse3, kInverse2);
        AVX2FF8T8VectorIFFTRadix4(
            value4, value5, value6, value7,
            kInverse5, kInverse7, kInverse6);
        AVX2FF8T8VectorIFFTDistance4(
            value0, value1, value2, value3,
            value4, value5, value6, value7, kInverse4);

        AVX2FF8T8VectorFFTRadix4Distance2(
            value0, value1, value2, value3,
            value4, value5, value6, value7,
            kForward2, kForward6, kForward4);
        AVX2FF8T8VectorFFT2(value0, value1, kForward1);
        AVX2FF8T8VectorFFT2(value2, value3, kForward3);
        AVX2FF8T8VectorFFT2(value4, value5, kForward5);
        AVX2FF8T8VectorFFT2(value6, value7, kForward7);

        AVX2FF8T8VectorMulAdd(value0, tail, kTail0);
        AVX2FF8T8VectorMulAdd(value1, tail, kTail1);
        AVX2FF8T8VectorMulAdd(value2, tail, kTail2);
        AVX2FF8T8VectorMulAdd(value3, tail, kTail3);
        AVX2FF8T8VectorMulAdd(value4, tail, kTail4);
        if (OutputCount >= 6)
            AVX2FF8T8VectorMulAdd(value5, tail, kTail5);
        if (OutputCount >= 7)
            AVX2FF8T8VectorMulAdd(value6, tail, kTail6);
        if (OutputCount >= 8)
            AVX2FF8T8VectorMulAdd(value7, tail, kTail7);

        _mm256_storeu_si256(
            reinterpret_cast<__m256i*>(output[0] + offset), value0);
        _mm256_storeu_si256(
            reinterpret_cast<__m256i*>(output[1] + offset), value1);
        _mm256_storeu_si256(
            reinterpret_cast<__m256i*>(output[2] + offset), value2);
        _mm256_storeu_si256(
            reinterpret_cast<__m256i*>(output[3] + offset), value3);
        _mm256_storeu_si256(
            reinterpret_cast<__m256i*>(output[4] + offset), value4);
        if (OutputCount >= 6)
            _mm256_storeu_si256(
                reinterpret_cast<__m256i*>(output[5] + offset), value5);
        if (OutputCount >= 7)
            _mm256_storeu_si256(
                reinterpret_cast<__m256i*>(output[6] + offset), value6);
        if (OutputCount >= 8)
            _mm256_storeu_si256(
                reinterpret_cast<__m256i*>(output[7] + offset), value7);
    }
}

/*
    Exact K=10..11/R=5..8 companions to the K=9 tail circuit.  The first
    eight sources are a complete T=8 message block.  Sources eight onward are
    parent coordinates 16 onward; their generator-column logarithms follow
    directly from

      L_x(y) = Z(y) / ((y + x) Z'(x)),
      Z(t) = product_{s=8}^{31} (t + s).

    Keeping the two tail loads after the first-block transform limits live
    vector pressure to nine values during each multiply-add sequence.
*/
template<unsigned TailCount, unsigned OutputCount>
static LEO2_AVX2_T8_ENTRY void AVX2FF8HighEncodeT8TailB256(
    const void* const* data,
    void* const* work)
{
    static_assert(TailCount >= 2 && TailCount <= 3,
        "T8 B256 tail count must be two or three");
    static_assert(OutputCount >= 5 && OutputCount <= 8,
        "K10 tail terminal output count must be between five and eight");
    static const uint16_t kInverse1 = 153;
    static const uint16_t kInverse2 = 17;
    static const uint16_t kInverse3 = 102;
    static const uint16_t kInverse4 = 85;
    static const uint16_t kInverse5 = 51;
    static const uint16_t kInverse6 = 34;
    static const uint16_t kInverse7 = 187;
    static const uint16_t kForward1 = 255;
    static const uint16_t kForward2 = 255;
    static const uint16_t kForward3 = 85;
    static const uint16_t kForward4 = 255;
    static const uint16_t kForward5 = 17;
    static const uint16_t kForward6 = 85;
    static const uint16_t kForward7 = 34;
    static const uint16_t kTail0[8] = {
        121, 151, 78, 228, 229, 94, 57, 147
    };
    static const uint16_t kTail1[8] = {
        151, 121, 228, 78, 94, 229, 147, 57
    };
    static const uint16_t kTail2[8] = {
        78, 228, 121, 151, 57, 147, 229, 94
    };

    /* Fixed extent keeps the compile-time-false TailCount=2 third-tail block
       syntactically in bounds under aggressive static analysis. */
    const uint8_t* input[11] = {};
    uint8_t* output[OutputCount];
    for (unsigned lane = 0; lane < 8 + TailCount; ++lane)
        input[lane] = static_cast<const uint8_t*>(data[lane]);
    for (unsigned lane = 0; lane < OutputCount; ++lane)
        output[lane] = static_cast<uint8_t*>(work[lane]);

    for (uint64_t offset = 0; offset < 256; offset += 32)
    {
        __m256i value0 = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(input[0] + offset));
        __m256i value1 = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(input[1] + offset));
        __m256i value2 = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(input[2] + offset));
        __m256i value3 = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(input[3] + offset));
        __m256i value4 = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(input[4] + offset));
        __m256i value5 = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(input[5] + offset));
        __m256i value6 = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(input[6] + offset));
        __m256i value7 = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(input[7] + offset));

        AVX2FF8T8VectorIFFTRadix4(
            value0, value1, value2, value3,
            kInverse1, kInverse3, kInverse2);
        AVX2FF8T8VectorIFFTRadix4(
            value4, value5, value6, value7,
            kInverse5, kInverse7, kInverse6);
        AVX2FF8T8VectorIFFTDistance4(
            value0, value1, value2, value3,
            value4, value5, value6, value7, kInverse4);

        AVX2FF8T8VectorFFTRadix4Distance2(
            value0, value1, value2, value3,
            value4, value5, value6, value7,
            kForward2, kForward6, kForward4);
        AVX2FF8T8VectorFFT2(value0, value1, kForward1);
        AVX2FF8T8VectorFFT2(value2, value3, kForward3);
        AVX2FF8T8VectorFFT2(value4, value5, kForward5);
        AVX2FF8T8VectorFFT2(value6, value7, kForward7);

        const __m256i tail0 = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(input[8] + offset));
        AVX2FF8T8VectorMulAdd(value0, tail0, kTail0[0]);
        AVX2FF8T8VectorMulAdd(value1, tail0, kTail0[1]);
        AVX2FF8T8VectorMulAdd(value2, tail0, kTail0[2]);
        AVX2FF8T8VectorMulAdd(value3, tail0, kTail0[3]);
        AVX2FF8T8VectorMulAdd(value4, tail0, kTail0[4]);
        if (OutputCount >= 6)
            AVX2FF8T8VectorMulAdd(value5, tail0, kTail0[5]);
        if (OutputCount >= 7)
            AVX2FF8T8VectorMulAdd(value6, tail0, kTail0[6]);
        if (OutputCount >= 8)
            AVX2FF8T8VectorMulAdd(value7, tail0, kTail0[7]);

        const __m256i tail1 = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(input[9] + offset));
        AVX2FF8T8VectorMulAdd(value0, tail1, kTail1[0]);
        AVX2FF8T8VectorMulAdd(value1, tail1, kTail1[1]);
        AVX2FF8T8VectorMulAdd(value2, tail1, kTail1[2]);
        AVX2FF8T8VectorMulAdd(value3, tail1, kTail1[3]);
        AVX2FF8T8VectorMulAdd(value4, tail1, kTail1[4]);
        if (OutputCount >= 6)
            AVX2FF8T8VectorMulAdd(value5, tail1, kTail1[5]);
        if (OutputCount >= 7)
            AVX2FF8T8VectorMulAdd(value6, tail1, kTail1[6]);
        if (OutputCount >= 8)
            AVX2FF8T8VectorMulAdd(value7, tail1, kTail1[7]);

        if (TailCount >= 3)
        {
            const __m256i tail2 = _mm256_loadu_si256(
                reinterpret_cast<const __m256i*>(input[10] + offset));
            AVX2FF8T8VectorMulAdd(value0, tail2, kTail2[0]);
            AVX2FF8T8VectorMulAdd(value1, tail2, kTail2[1]);
            AVX2FF8T8VectorMulAdd(value2, tail2, kTail2[2]);
            AVX2FF8T8VectorMulAdd(value3, tail2, kTail2[3]);
            AVX2FF8T8VectorMulAdd(value4, tail2, kTail2[4]);
            if (OutputCount >= 6)
                AVX2FF8T8VectorMulAdd(value5, tail2, kTail2[5]);
            if (OutputCount >= 7)
                AVX2FF8T8VectorMulAdd(value6, tail2, kTail2[6]);
            if (OutputCount >= 8)
                AVX2FF8T8VectorMulAdd(value7, tail2, kTail2[7]);
        }

        _mm256_storeu_si256(
            reinterpret_cast<__m256i*>(output[0] + offset), value0);
        _mm256_storeu_si256(
            reinterpret_cast<__m256i*>(output[1] + offset), value1);
        _mm256_storeu_si256(
            reinterpret_cast<__m256i*>(output[2] + offset), value2);
        _mm256_storeu_si256(
            reinterpret_cast<__m256i*>(output[3] + offset), value3);
        _mm256_storeu_si256(
            reinterpret_cast<__m256i*>(output[4] + offset), value4);
        if (OutputCount >= 6)
            _mm256_storeu_si256(
                reinterpret_cast<__m256i*>(output[5] + offset), value5);
        if (OutputCount >= 7)
            _mm256_storeu_si256(
                reinterpret_cast<__m256i*>(output[6] + offset), value6);
        if (OutputCount >= 8)
            _mm256_storeu_si256(
                reinterpret_cast<__m256i*>(output[7] + offset), value7);
    }
}

/*
    R=8 comparison form: keep each pair's four nibble tables live for the
    complete shard instead of reloading them for every 32-byte transform
    slice.  The extra parity read/write pass is measured against the fused
    form before this path is retained.
*/
static LEO2_AVX2_T8_ENTRY void AVX2FF8HighEncodeK10R8T8B256Pairwise(
    const void* const* data,
    void* const* work)
{
    static const uint8_t kInverseSkew[8] = {
        255, 153, 17, 102, 85, 51, 34, 187
    };
    static const uint8_t kForwardSkew[8] = {
        0, 255, 255, 85, 255, 17, 85, 34
    };
    static const uint16_t kTail0[8] = {
        121, 151, 78, 228, 229, 94, 57, 147
    };
    static const uint16_t kTail1[8] = {
        151, 121, 228, 78, 94, 229, 147, 57
    };
    AVX2FF8HighEncodeT8Vector(
        data, work, false, kInverseSkew, kForwardSkew, 256);
    for (unsigned output = 0; output < 8; output += 2)
    {
        void* destinations[2] = { work[output], work[output + 1] };
        AVX2FF8MultiplyAdd2Sources2Outputs(
            destinations, data[8], data[9],
            kTail0 + output, kTail1 + output, 256);
    }
}
#endif // LEO2_EXPERIMENT_HIGH_T8_TWO_BLOCK_BINDING

/*
    Exact shortened/punctured T=8 tile for the K=5,R=5 boundary.  Unlike the
    padded callback, this keeps the three known-zero input lanes in registers
    and never writes the three punctured parity lanes.  The arithmetic is the
    same complete parent transform; only provably dead loads and stores are
    removed.
*/
static LEO2_AVX2_T8_ENTRY void AVX2FF8HighEncodeK5R5T8Multiple32(
    const void* const* data,
    void* const* work,
    uint64_t byte_count)
{
    /*
        These are the deterministic legacy-GF8 T=8 message/parity coset
        factors produced by FFTInitialize.  Making them integral constants
        lets the compiler remove zero-skew branches and address the fixed
        nibble tables directly in this versioned wire-profile kernel.
        Backend qualification independently evaluates the same factors through
        the generic reference transform.
    */
    static const uint16_t kInverse1 = 153;
    static const uint16_t kInverse2 = 17;
    static const uint16_t kInverse3 = 102;
    static const uint16_t kInverse4 = 85;
    static const uint16_t kInverse5 = 51;
    static const uint16_t kInverse6 = 34;
    static const uint16_t kForward5 = 17;
    LEO_DEBUG_ASSERT(byte_count != 0 && (byte_count & 31) == 0);

    const uint8_t* input0 = static_cast<const uint8_t*>(data[0]);
    const uint8_t* input1 = static_cast<const uint8_t*>(data[1]);
    const uint8_t* input2 = static_cast<const uint8_t*>(data[2]);
    const uint8_t* input3 = static_cast<const uint8_t*>(data[3]);
    const uint8_t* input4 = static_cast<const uint8_t*>(data[4]);
    uint8_t* output0 = static_cast<uint8_t*>(work[0]);
    uint8_t* output1 = static_cast<uint8_t*>(work[1]);
    uint8_t* output2 = static_cast<uint8_t*>(work[2]);
    uint8_t* output3 = static_cast<uint8_t*>(work[3]);
    uint8_t* output4 = static_cast<uint8_t*>(work[4]);

    for (uint64_t offset = 0; offset < byte_count; offset += 32)
    {
        __m256i value0 = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(input0 + offset));
        __m256i value1 = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(input1 + offset));
        __m256i value2 = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(input2 + offset));
        __m256i value3 = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(input3 + offset));
        __m256i value4 = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(input4 + offset));

        AVX2FF8T8VectorIFFTRadix4(
            value0, value1, value2, value3,
            kInverse1, kInverse3, kInverse2);

        /*
            The shortened upper inverse radix-4 maps [E,0,0,0] to

              [M119(E), M136(E), M238(E), E].

            Generate the same values with its two useful factors, then fuse
            each distance-4 inverse pair with the cancelling forward XOR:

              t = u ^ s; z = u ^ M85(t); b = t ^ z.

            The five requested outputs need only b0, b1 and their products
            with b2, b3.  This avoids materializing the three punctured
            forward outputs and, importantly, keeps the log-17 tables out of
            the eight-value live range.
        */
        __m256i value5 = value4;
        AVX2FF8T8VectorMulAdd(value5, value4, kInverse5);
        __m256i value6 = value5;
        __m256i value7 = value4;
        {
            const FF8NibbleTable& table = FF8Tables[kInverse6];
            const __m256i low_table = BroadcastTable(table.low);
            const __m256i high_table = BroadcastTable(table.high);
            AVX2FF8T8VectorMulAddPrepared(
                value6, value5, low_table, high_table);
            AVX2FF8T8VectorMulAddPrepared(
                value7, value4, low_table, high_table);
        }
        {
            const FF8NibbleTable& table = FF8Tables[kInverse4];
            const __m256i low_table = BroadcastTable(table.low);
            const __m256i high_table = BroadcastTable(table.high);

            AVX2FF8T8VectorXor(value6, value0);
            AVX2FF8T8VectorMulAddPrepared(
                value0, value6, low_table, high_table);
            AVX2FF8T8VectorXor(value6, value0);

            AVX2FF8T8VectorXor(value5, value2);
            AVX2FF8T8VectorMulAddPrepared(
                value2, value5, low_table, high_table);
            AVX2FF8T8VectorXor(value5, value2);

            AVX2FF8T8VectorXor(value7, value1);
            AVX2FF8T8VectorMulAddPrepared(
                value1, value7, low_table, high_table);
            AVX2FF8T8VectorXor(value7, value1);

            AVX2FF8T8VectorXor(value4, value3);
            AVX2FF8T8VectorMulAddPrepared(
                value3, value4, low_table, high_table);
            AVX2FF8T8VectorXor(value4, value3);

            AVX2FF8T8VectorXor(value2, value0);
            AVX2FF8T8VectorXor(value3, value1);
            AVX2FF8T8VectorMulAddPrepared(
                value2, value3, low_table, high_table);
            AVX2FF8T8VectorXor(value3, value2);
            AVX2FF8T8VectorXor(value1, value0);

            /*
                Scatter the lower four outputs before reducing the upper
                pair.  This deliberately ends their live ranges while the
                log-85 tables are still resident, leaving enough registers
                for the two upper reductions without a stack spill.
            */
            _mm256_storeu_si256(
                reinterpret_cast<__m256i*>(output0 + offset), value0);
            _mm256_storeu_si256(
                reinterpret_cast<__m256i*>(output1 + offset), value1);
            _mm256_storeu_si256(
                reinterpret_cast<__m256i*>(output2 + offset), value2);
            _mm256_storeu_si256(
                reinterpret_cast<__m256i*>(output3 + offset), value3);

            AVX2FF8T8VectorMulAddPrepared(
                value6, value5, low_table, high_table);
            AVX2FF8T8VectorMulAddPrepared(
                value7, value4, low_table, high_table);
        }
        AVX2FF8T8VectorMulAdd(value6, value7, kForward5);

        _mm256_storeu_si256(
            reinterpret_cast<__m256i*>(output4 + offset), value6);
    }
}

static LEO2_AVX2_T8_ENTRY void AVX2FF8HighEncodeK5R5T8Tail(
    const void* const* data,
    void* const* work,
    uint64_t byte_count)
{
    static const uint16_t kInverse1 = 153;
    static const uint16_t kInverse2 = 17;
    static const uint16_t kInverse3 = 102;
    static const uint16_t kInverse4 = 85;
    static const uint16_t kInverse5 = 51;
    static const uint16_t kInverse6 = 34;
    static const uint16_t kInverse7 = 187;
    static const uint16_t kForward1 = 255;
    static const uint16_t kForward2 = 255;
    static const uint16_t kForward3 = 85;
    static const uint16_t kForward4 = 255;
    static const uint16_t kForward5 = 17;
    static const uint16_t kForward6 = 85;

    const uint64_t vector_bytes = byte_count & ~UINT64_C(31);
    if (vector_bytes != 0)
        AVX2FF8HighEncodeK5R5T8Multiple32(data, work, vector_bytes);

    const uint8_t* input[5];
    uint8_t* output[5];
    for (unsigned lane = 0; lane < 5; ++lane)
    {
        input[lane] = static_cast<const uint8_t*>(data[lane]);
        output[lane] = static_cast<uint8_t*>(work[lane]);
    }

    for (uint64_t offset = vector_bytes; offset < byte_count; ++offset)
    {
        uint8_t values[8] = {};
        for (unsigned lane = 0; lane < 5; ++lane)
            values[lane] = input[lane][offset];

        AVX2FF8T8VectorScalarIFFT2(
            values[0], values[1], kInverse1);
        AVX2FF8T8VectorScalarIFFT2(
            values[2], values[3], kInverse3);
        AVX2FF8T8VectorScalarIFFT2(
            values[0], values[2], kInverse2);
        AVX2FF8T8VectorScalarIFFT2(
            values[1], values[3], kInverse2);
        AVX2FF8T8VectorScalarIFFT2(
            values[4], values[5], kInverse5);
        AVX2FF8T8VectorScalarIFFT2(
            values[6], values[7], kInverse7);
        AVX2FF8T8VectorScalarIFFT2(
            values[4], values[6], kInverse6);
        AVX2FF8T8VectorScalarIFFT2(
            values[5], values[7], kInverse6);
        for (unsigned lane = 0; lane < 4; ++lane)
            AVX2FF8T8VectorScalarIFFT2(
                values[lane], values[lane + 4], kInverse4);

        AVX2FF8T8VectorScalarFFT2(
            values[0], values[4], kForward4);
        AVX2FF8T8VectorScalarFFT2(
            values[1], values[5], kForward4);
        AVX2FF8T8VectorScalarFFT2(
            values[2], values[6], kForward4);
        AVX2FF8T8VectorScalarFFT2(
            values[3], values[7], kForward4);
        AVX2FF8T8VectorScalarFFT2(
            values[0], values[2], kForward2);
        AVX2FF8T8VectorScalarFFT2(
            values[1], values[3], kForward2);
        AVX2FF8T8VectorScalarFFT2(
            values[4], values[6], kForward6);
        AVX2FF8T8VectorScalarFFT2(
            values[5], values[7], kForward6);
        AVX2FF8T8VectorScalarFFT2(
            values[0], values[1], kForward1);
        AVX2FF8T8VectorScalarFFT2(
            values[2], values[3], kForward3);
        AVX2FF8T8VectorScalarFFT2(
            values[4], values[5], kForward5);

        for (unsigned lane = 0; lane < 5; ++lane)
            output[lane][offset] = values[lane];
    }
}

static LEO2_AVX2_T8_ENTRY void AVX2FF8HighEncodeK5R5T8Vector(
    const void* const* data,
    void* const* work,
    const uint8_t* inverse_skew,
    const uint8_t* forward_skew,
    uint64_t byte_count)
{
    LEO_DEBUG_ASSERT(
        inverse_skew[0] == 255 &&
        inverse_skew[1] == 153 &&
        inverse_skew[2] == 17 &&
        inverse_skew[3] == 102 &&
        inverse_skew[4] == 85 &&
        inverse_skew[5] == 51 &&
        inverse_skew[6] == 34 &&
        inverse_skew[7] == 187);
    LEO_DEBUG_ASSERT(
        forward_skew[0] == 0 &&
        forward_skew[1] == 255 &&
        forward_skew[2] == 255 &&
        forward_skew[3] == 85 &&
        forward_skew[4] == 255 &&
        forward_skew[5] == 17 &&
        forward_skew[6] == 85 &&
        forward_skew[7] == 34);
    (void)inverse_skew;
    (void)forward_skew;

    if (byte_count == 0)
        return;
    if ((byte_count & 31) == 0)
    {
        AVX2FF8HighEncodeK5R5T8Multiple32(data, work, byte_count);
        return;
    }
    AVX2FF8HighEncodeK5R5T8Tail(data, work, byte_count);
}

/*
    Exact shortened/punctured companions for K=6 through K=8.  The complete T=8
    parent transform is unchanged, but known-zero source suffixes are formed
    in registers and punctured parity suffixes are never stored.  OutputCount
    permits the K=6/R=5 prefix; dead final-forward work is removed when this
    template is instantiated.
*/
template<unsigned ActiveCount, unsigned OutputCount = ActiveCount>
static LEO2_AVX2_T8_ENTRY void AVX2FF8HighEncodePartialT8Vector(
    const void* const* data,
    void* const* work,
    const uint8_t* inverse_skew,
    const uint8_t* forward_skew,
    uint64_t byte_count)
{
    static_assert(ActiveCount >= 6 && ActiveCount <= 8,
        "exact T=8 partial kernel supports K=6 through K=8");
    static_assert(OutputCount >= 5 && OutputCount <= ActiveCount,
        "exact T=8 partial output prefix is invalid");
    static const uint16_t kInverse1 = 153;
    static const uint16_t kInverse2 = 17;
    static const uint16_t kInverse3 = 102;
    static const uint16_t kInverse4 = 85;
    static const uint16_t kInverse5 = 51;
    static const uint16_t kInverse6 = 34;
    static const uint16_t kInverse7 = 187;
    static const uint16_t kForward1 = 255;
    static const uint16_t kForward2 = 255;
    static const uint16_t kForward3 = 85;
    static const uint16_t kForward4 = 255;
    static const uint16_t kForward5 = 17;
    static const uint16_t kForward6 = 85;
    static const uint16_t kForward7 = 34;
    LEO_DEBUG_ASSERT(
        inverse_skew[0] == 255 &&
        inverse_skew[1] == kInverse1 &&
        inverse_skew[2] == kInverse2 &&
        inverse_skew[3] == kInverse3 &&
        inverse_skew[4] == kInverse4 &&
        inverse_skew[5] == kInverse5 &&
        inverse_skew[6] == kInverse6 &&
        inverse_skew[7] == kInverse7);
    LEO_DEBUG_ASSERT(
        forward_skew[0] == 0 &&
        forward_skew[1] == kForward1 &&
        forward_skew[2] == kForward2 &&
        forward_skew[3] == kForward3 &&
        forward_skew[4] == kForward4 &&
        forward_skew[5] == kForward5 &&
        forward_skew[6] == kForward6 &&
        forward_skew[7] == kForward7);
    (void)inverse_skew;
    (void)forward_skew;

    const uint8_t* input[ActiveCount];
    for (unsigned lane = 0; lane < ActiveCount; ++lane)
        input[lane] = static_cast<const uint8_t*>(data[lane]);
    uint8_t* output0 = static_cast<uint8_t*>(work[0]);
    uint8_t* output1 = static_cast<uint8_t*>(work[1]);
    uint8_t* output2 = static_cast<uint8_t*>(work[2]);
    uint8_t* output3 = static_cast<uint8_t*>(work[3]);
    uint8_t* output4 = static_cast<uint8_t*>(work[4]);
    uint8_t* output5 = OutputCount >= 6
        ? static_cast<uint8_t*>(work[5]) : NULL;
    uint8_t* output6 = OutputCount >= 7
        ? static_cast<uint8_t*>(work[6]) : NULL;
    uint8_t* output7 = OutputCount >= 8
        ? static_cast<uint8_t*>(work[7]) : NULL;

    const uint64_t vector_bytes = byte_count & ~UINT64_C(31);
    for (uint64_t offset = 0; offset < vector_bytes; offset += 32)
    {
        __m256i value0 = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(input[0] + offset));
        __m256i value1 = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(input[1] + offset));
        __m256i value2 = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(input[2] + offset));
        __m256i value3 = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(input[3] + offset));
        __m256i value4 = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(input[4] + offset));
        __m256i value5 = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(input[5] + offset));
        __m256i value6 = _mm256_setzero_si256();
        if (ActiveCount >= 7)
            value6 = _mm256_loadu_si256(
                reinterpret_cast<const __m256i*>(input[6] + offset));
        __m256i value7 = _mm256_setzero_si256();
        if (ActiveCount == 8)
            value7 = _mm256_loadu_si256(
                reinterpret_cast<const __m256i*>(input[7] + offset));

        AVX2FF8T8VectorIFFTRadix4(
            value0, value1, value2, value3,
            kInverse1, kInverse3, kInverse2);
        AVX2FF8T8VectorIFFTRadix4(
            value4, value5, value6, value7,
            kInverse5, kInverse7, kInverse6);
        AVX2FF8T8VectorIFFTDistance4(
            value0, value1, value2, value3,
            value4, value5, value6, value7, kInverse4);

        AVX2FF8T8VectorFFTRadix4Distance2(
            value0, value1, value2, value3,
            value4, value5, value6, value7,
            kForward2, kForward6, kForward4);
        AVX2FF8T8VectorFFT2(value0, value1, kForward1);
        AVX2FF8T8VectorFFT2(value2, value3, kForward3);
        AVX2FF8T8VectorFFT2(value4, value5, kForward5);
        AVX2FF8T8VectorFFT2(value6, value7, kForward7);

        _mm256_storeu_si256(
            reinterpret_cast<__m256i*>(output0 + offset), value0);
        _mm256_storeu_si256(
            reinterpret_cast<__m256i*>(output1 + offset), value1);
        _mm256_storeu_si256(
            reinterpret_cast<__m256i*>(output2 + offset), value2);
        _mm256_storeu_si256(
            reinterpret_cast<__m256i*>(output3 + offset), value3);
        _mm256_storeu_si256(
            reinterpret_cast<__m256i*>(output4 + offset), value4);
        if (OutputCount >= 6)
            _mm256_storeu_si256(
                reinterpret_cast<__m256i*>(output5 + offset), value5);
        if (OutputCount >= 7)
            _mm256_storeu_si256(
                reinterpret_cast<__m256i*>(output6 + offset), value6);
        if (OutputCount >= 8)
            _mm256_storeu_si256(
                reinterpret_cast<__m256i*>(output7 + offset), value7);
    }

    for (uint64_t offset = vector_bytes; offset < byte_count; ++offset)
    {
        uint8_t values[8] = {};
        for (unsigned lane = 0; lane < ActiveCount; ++lane)
            values[lane] = input[lane][offset];

        AVX2FF8T8VectorScalarIFFT2(
            values[0], values[1], kInverse1);
        AVX2FF8T8VectorScalarIFFT2(
            values[2], values[3], kInverse3);
        AVX2FF8T8VectorScalarIFFT2(
            values[0], values[2], kInverse2);
        AVX2FF8T8VectorScalarIFFT2(
            values[1], values[3], kInverse2);
        AVX2FF8T8VectorScalarIFFT2(
            values[4], values[5], kInverse5);
        AVX2FF8T8VectorScalarIFFT2(
            values[6], values[7], kInverse7);
        AVX2FF8T8VectorScalarIFFT2(
            values[4], values[6], kInverse6);
        AVX2FF8T8VectorScalarIFFT2(
            values[5], values[7], kInverse6);
        for (unsigned lane = 0; lane < 4; ++lane)
            AVX2FF8T8VectorScalarIFFT2(
                values[lane], values[lane + 4], kInverse4);

        AVX2FF8T8VectorScalarFFT2(
            values[0], values[4], kForward4);
        AVX2FF8T8VectorScalarFFT2(
            values[1], values[5], kForward4);
        AVX2FF8T8VectorScalarFFT2(
            values[2], values[6], kForward4);
        AVX2FF8T8VectorScalarFFT2(
            values[3], values[7], kForward4);
        AVX2FF8T8VectorScalarFFT2(
            values[0], values[2], kForward2);
        AVX2FF8T8VectorScalarFFT2(
            values[1], values[3], kForward2);
        AVX2FF8T8VectorScalarFFT2(
            values[4], values[6], kForward6);
        AVX2FF8T8VectorScalarFFT2(
            values[5], values[7], kForward6);
        AVX2FF8T8VectorScalarFFT2(
            values[0], values[1], kForward1);
        AVX2FF8T8VectorScalarFFT2(
            values[2], values[3], kForward3);
        AVX2FF8T8VectorScalarFFT2(
            values[4], values[5], kForward5);
        AVX2FF8T8VectorScalarFFT2(
            values[6], values[7], kForward7);

        for (unsigned lane = 0; lane < OutputCount; ++lane)
            static_cast<uint8_t*>(work[lane])[offset] = values[lane];
    }
}

static LEO2_AVX2_T8_ENTRY void AVX2FF8HighEncodeOneBlockT8Vector(
    const void* const* data,
    void* const* work,
    uint32_t side_and_flags,
    const uint8_t* inverse_skew,
    const uint8_t* forward_skew,
    uint64_t byte_count)
{
    const bool shortened =
        (side_and_flags & kFF8HighEncodeShortenedInput) != 0;
    const bool k5_r5_partial =
        (side_and_flags & kFF8HighEncodeK5R5Partial) != 0;
    const bool k6_r6_partial =
        (side_and_flags & kFF8HighEncodeK6R6Partial) != 0;
    const bool k7_r7_partial =
        (side_and_flags & kFF8HighEncodeK7R7Partial) != 0;
    const bool k6_r5_partial =
        (side_and_flags & kFF8HighEncodeK6R5Partial) != 0;
    const bool k7_r5_partial =
        (side_and_flags & kFF8HighEncodeK7R5Partial) != 0;
    const bool k8_r5_partial =
        (side_and_flags & kFF8HighEncodeK8R5Partial) != 0;
#if LEO2_EXPERIMENT_HIGH_T8_TWO_BLOCK_BINDING
    const bool k9_tail =
        (side_and_flags & kFF8HighEncodeK9Tail) != 0;
    const uint32_t k9_output_count =
        (side_and_flags & kFF8HighEncodeK9OutputCountMask) >>
            kFF8HighEncodeK9OutputCountShift;
    const uint32_t side = side_and_flags &
        ~(kFF8HighEncodeShortenedInput |
          kFF8HighEncodeK5R5Partial |
          kFF8HighEncodeK6R6Partial |
          kFF8HighEncodeK7R7Partial |
          kFF8HighEncodeK6R5Partial |
          kFF8HighEncodeK7R5Partial |
          kFF8HighEncodeK8R5Partial |
          kFF8HighEncodeK9Tail |
          kFF8HighEncodeK9OutputCountMask);
#else
    const bool k9_tail = false;
    const uint32_t k9_output_count = 0;
    const uint32_t side =
        side_and_flags &
        ~(kFF8HighEncodeShortenedInput |
          kFF8HighEncodeK5R5Partial |
          kFF8HighEncodeK6R6Partial |
          kFF8HighEncodeK7R7Partial |
          kFF8HighEncodeK6R5Partial |
          kFF8HighEncodeK7R5Partial |
          kFF8HighEncodeK8R5Partial);
#endif
    LEO_DEBUG_ASSERT(side == 8);
    if (side != 8)
        return;
    const bool k9_output_count_valid =
        k9_output_count >= 5 && k9_output_count <= 8;
    const unsigned partial_count =
        static_cast<unsigned>(k5_r5_partial) +
        static_cast<unsigned>(k6_r6_partial) +
        static_cast<unsigned>(k7_r7_partial) +
        static_cast<unsigned>(k6_r5_partial) +
        static_cast<unsigned>(k7_r5_partial) +
        static_cast<unsigned>(k8_r5_partial);
    if ((shortened && (partial_count != 0 || k9_tail)) ||
        partial_count > 1 || (partial_count != 0 && k9_tail) ||
        (k9_tail ? !k9_output_count_valid : k9_output_count != 0))
        return;
    LEO_DEBUG_ASSERT(
        static_cast<unsigned>(shortened) +
        partial_count +
        static_cast<unsigned>(k9_tail) <= 1U);
#if LEO2_EXPERIMENT_HIGH_T8_TWO_BLOCK_BINDING
    if (k9_tail)
    {
        LEO_DEBUG_ASSERT(
            inverse_skew[0] == 255 &&
            inverse_skew[1] == 153 &&
            inverse_skew[2] == 17 &&
            inverse_skew[3] == 102 &&
            inverse_skew[4] == 85 &&
            inverse_skew[5] == 51 &&
            inverse_skew[6] == 34 &&
            inverse_skew[7] == 187);
        LEO_DEBUG_ASSERT(
            forward_skew[0] == 0 &&
            forward_skew[1] == 255 &&
            forward_skew[2] == 255 &&
            forward_skew[3] == 85 &&
            forward_skew[4] == 255 &&
            forward_skew[5] == 17 &&
            forward_skew[6] == 85 &&
            forward_skew[7] == 34);
        (void)inverse_skew;
        (void)forward_skew;
        switch (k9_output_count)
        {
        case 5:
            AVX2FF8HighEncodeK9R5T8Multiple32(data, work, byte_count);
            break;
        case 6:
            AVX2FF8HighEncodeK9T8Multiple32<6>(data, work, byte_count);
            break;
        case 7:
            AVX2FF8HighEncodeK9T8Multiple32<7>(data, work, byte_count);
            break;
        case 8:
            AVX2FF8HighEncodeK9T8Multiple32<8>(data, work, byte_count);
            break;
        default:
            LEO_DEBUG_ASSERT(false);
            break;
        }
        return;
    }
#endif
    if (k5_r5_partial)
    {
        LEO_DEBUG_ASSERT(!shortened);
        AVX2FF8HighEncodeK5R5T8Vector(
            data, work, inverse_skew, forward_skew, byte_count);
        return;
    }
    if (k6_r6_partial)
    {
        AVX2FF8HighEncodePartialT8Vector<6>(
            data, work, inverse_skew, forward_skew, byte_count);
        return;
    }
    if (k6_r5_partial)
    {
        AVX2FF8HighEncodePartialT8Vector<6, 5>(
            data, work, inverse_skew, forward_skew, byte_count);
        return;
    }
    if (k7_r7_partial)
    {
        AVX2FF8HighEncodePartialT8Vector<7>(
            data, work, inverse_skew, forward_skew, byte_count);
        return;
    }
    if (k7_r5_partial)
    {
        AVX2FF8HighEncodePartialT8Vector<7, 5>(
            data, work, inverse_skew, forward_skew, byte_count);
        return;
    }
    if (k8_r5_partial)
    {
        AVX2FF8HighEncodePartialT8Vector<8, 5>(
            data, work, inverse_skew, forward_skew, byte_count);
        return;
    }
    AVX2FF8HighEncodeT8Vector(
        data, work, shortened, inverse_skew, forward_skew, byte_count);
}

#if LEO2_EXPERIMENT_HIGH_T8_TWO_BLOCK_BINDING
static LEO_FORCE_INLINE void AVX2FF8T8InverseValues(
    __m256i& value0,
    __m256i& value1,
    __m256i& value2,
    __m256i& value3,
    __m256i& value4,
    __m256i& value5,
    __m256i& value6,
    __m256i& value7,
    const uint8_t* inverse_skew)
{
    AVX2FF8T8VectorIFFTRadix4(
        value0, value1, value2, value3,
        inverse_skew[1], inverse_skew[3], inverse_skew[2]);
    AVX2FF8T8VectorIFFTRadix4(
        value4, value5, value6, value7,
        inverse_skew[5], inverse_skew[7], inverse_skew[6]);
    AVX2FF8T8VectorIFFTDistance4(
        value0, value1, value2, value3,
        value4, value5, value6, value7, inverse_skew[4]);
}

static LEO_FORCE_INLINE void AVX2FF8T8ForwardValues(
    __m256i& value0,
    __m256i& value1,
    __m256i& value2,
    __m256i& value3,
    __m256i& value4,
    __m256i& value5,
    __m256i& value6,
    __m256i& value7,
    const uint8_t* forward_skew)
{
    AVX2FF8T8VectorFFTRadix4Distance2(
        value0, value1, value2, value3,
        value4, value5, value6, value7,
        forward_skew[2], forward_skew[6], forward_skew[4]);
    AVX2FF8T8VectorFFT2(value0, value1, forward_skew[1]);
    AVX2FF8T8VectorFFT2(value2, value3, forward_skew[3]);
    AVX2FF8T8VectorFFT2(value4, value5, forward_skew[5]);
    AVX2FF8T8VectorFFT2(value6, value7, forward_skew[7]);
}

static LEO_FORCE_INLINE void AVX2FF8HighEncodeTwoBlocksT8Vectors(
    const void* const* data,
    void* const* work,
    const uint8_t* first_inverse_skew,
    const uint8_t* second_inverse_skew,
    const uint8_t* forward_skew,
    uint64_t first_offset,
    uint64_t byte_count)
{
    LEO_DEBUG_ASSERT(
        byte_count >= 32 && byte_count <= 1024 &&
        (byte_count & 31U) == 0);
    if (byte_count < 32 || byte_count > 1024 ||
        (byte_count & 31U) != 0)
        return;

    const uint64_t end_offset = first_offset + byte_count;
    for (uint64_t offset = first_offset;
         offset < end_offset; offset += 32)
    {
        __m256i value0 = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(
                static_cast<const uint8_t*>(data[0]) + offset));
        __m256i value1 = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(
                static_cast<const uint8_t*>(data[1]) + offset));
        __m256i value2 = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(
                static_cast<const uint8_t*>(data[2]) + offset));
        __m256i value3 = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(
                static_cast<const uint8_t*>(data[3]) + offset));
        __m256i value4 = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(
                static_cast<const uint8_t*>(data[4]) + offset));
        __m256i value5 = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(
                static_cast<const uint8_t*>(data[5]) + offset));
        __m256i value6 = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(
                static_cast<const uint8_t*>(data[6]) + offset));
        __m256i value7 = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(
                static_cast<const uint8_t*>(data[7]) + offset));
        AVX2FF8T8InverseValues(
            value0, value1, value2, value3,
            value4, value5, value6, value7, first_inverse_skew);
        _mm256_storeu_si256(
            reinterpret_cast<__m256i*>(
                static_cast<uint8_t*>(work[0]) + offset), value0);
        _mm256_storeu_si256(
            reinterpret_cast<__m256i*>(
                static_cast<uint8_t*>(work[1]) + offset), value1);
        _mm256_storeu_si256(
            reinterpret_cast<__m256i*>(
                static_cast<uint8_t*>(work[2]) + offset), value2);
        _mm256_storeu_si256(
            reinterpret_cast<__m256i*>(
                static_cast<uint8_t*>(work[3]) + offset), value3);
        _mm256_storeu_si256(
            reinterpret_cast<__m256i*>(
                static_cast<uint8_t*>(work[4]) + offset), value4);
        _mm256_storeu_si256(
            reinterpret_cast<__m256i*>(
                static_cast<uint8_t*>(work[5]) + offset), value5);
        _mm256_storeu_si256(
            reinterpret_cast<__m256i*>(
                static_cast<uint8_t*>(work[6]) + offset), value6);
        _mm256_storeu_si256(
            reinterpret_cast<__m256i*>(
                static_cast<uint8_t*>(work[7]) + offset), value7);

        value0 = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(
                static_cast<const uint8_t*>(data[8]) + offset));
        value1 = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(
                static_cast<const uint8_t*>(data[9]) + offset));
        value2 = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(
                static_cast<const uint8_t*>(data[10]) + offset));
        value3 = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(
                static_cast<const uint8_t*>(data[11]) + offset));
        value4 = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(
                static_cast<const uint8_t*>(data[12]) + offset));
        value5 = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(
                static_cast<const uint8_t*>(data[13]) + offset));
        value6 = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(
                static_cast<const uint8_t*>(data[14]) + offset));
        value7 = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(
                static_cast<const uint8_t*>(data[15]) + offset));
        AVX2FF8T8InverseValues(
            value0, value1, value2, value3,
            value4, value5, value6, value7, second_inverse_skew);
        value0 = _mm256_xor_si256(value0, _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(
                static_cast<const uint8_t*>(work[0]) + offset)));
        value1 = _mm256_xor_si256(value1, _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(
                static_cast<const uint8_t*>(work[1]) + offset)));
        value2 = _mm256_xor_si256(value2, _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(
                static_cast<const uint8_t*>(work[2]) + offset)));
        value3 = _mm256_xor_si256(value3, _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(
                static_cast<const uint8_t*>(work[3]) + offset)));
        value4 = _mm256_xor_si256(value4, _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(
                static_cast<const uint8_t*>(work[4]) + offset)));
        value5 = _mm256_xor_si256(value5, _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(
                static_cast<const uint8_t*>(work[5]) + offset)));
        value6 = _mm256_xor_si256(value6, _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(
                static_cast<const uint8_t*>(work[6]) + offset)));
        value7 = _mm256_xor_si256(value7, _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(
                static_cast<const uint8_t*>(work[7]) + offset)));
        AVX2FF8T8ForwardValues(
            value0, value1, value2, value3,
            value4, value5, value6, value7, forward_skew);
        _mm256_storeu_si256(
            reinterpret_cast<__m256i*>(
                static_cast<uint8_t*>(work[0]) + offset), value0);
        _mm256_storeu_si256(
            reinterpret_cast<__m256i*>(
                static_cast<uint8_t*>(work[1]) + offset), value1);
        _mm256_storeu_si256(
            reinterpret_cast<__m256i*>(
                static_cast<uint8_t*>(work[2]) + offset), value2);
        _mm256_storeu_si256(
            reinterpret_cast<__m256i*>(
                static_cast<uint8_t*>(work[3]) + offset), value3);
        _mm256_storeu_si256(
            reinterpret_cast<__m256i*>(
                static_cast<uint8_t*>(work[4]) + offset), value4);
        _mm256_storeu_si256(
            reinterpret_cast<__m256i*>(
                static_cast<uint8_t*>(work[5]) + offset), value5);
        _mm256_storeu_si256(
            reinterpret_cast<__m256i*>(
                static_cast<uint8_t*>(work[6]) + offset), value6);
        _mm256_storeu_si256(
            reinterpret_cast<__m256i*>(
                static_cast<uint8_t*>(work[7]) + offset), value7);
    }
}

/*
    Keep one copy of the exact-vector arithmetic for the tiny entry point.
    The mature callback still inlines its loop exactly as before; the bounded
    tiny wrapper pays at most two calls to avoid triplicating a large nibble-
    table transform in its instruction footprint.
*/
static LEO2_AVX2_T8_ENTRY void AVX2FF8HighEncodeTwoBlocksT8TinyVector32(
    const void* const* data,
    void* const* work,
    const uint8_t* first_inverse_skew,
    const uint8_t* second_inverse_skew,
    const uint8_t* forward_skew,
    uint64_t offset)
{
    AVX2FF8HighEncodeTwoBlocksT8Vectors(
        data, work, first_inverse_skew, second_inverse_skew,
        forward_skew, offset, 32);
}

static LEO2_AVX2_T8_ENTRY void AVX2FF8HighEncodeTwoBlocksT8(
    const void* const* data,
    void* const* work,
    const uint8_t* first_inverse_skew,
    const uint8_t* second_inverse_skew,
    const uint8_t* forward_skew,
    uint64_t byte_count)
{
    LEO_DEBUG_ASSERT(byte_count >= 64 && byte_count <= 1024);
    if (byte_count < 64 || byte_count > 1024)
        return;

    const uint64_t complete_bytes = byte_count & ~UINT64_C(31);
    AVX2FF8HighEncodeTwoBlocksT8Vectors(
        data, work, first_inverse_skew, second_inverse_skew,
        forward_skew, 0, complete_bytes);
    if (complete_bytes != byte_count)
    {
        /*
            Recompute the final in-range vector.  It may overlap the preceding
            vector, but every byte lane is independent and therefore produces
            the same output.  This avoids any out-of-range load or tail copy.
        */
        AVX2FF8HighEncodeTwoBlocksT8TinyVector32(
            data, work, first_inverse_skew, second_inverse_skew,
            forward_skew, byte_count - 32);
    }
}

static LEO2_AVX2_T8_ENTRY void AVX2FF8HighEncodeTwoBlocksT8Tiny(
    const void* const* data,
    void* const* work,
    const uint8_t* first_inverse_skew,
    const uint8_t* second_inverse_skew,
    const uint8_t* forward_skew,
    uint64_t byte_count)
{
    LEO_DEBUG_ASSERT(byte_count >= 1 && byte_count < 64);
    if (byte_count == 0 || byte_count >= 64)
        return;

    if (byte_count >= 32)
    {
        AVX2FF8HighEncodeTwoBlocksT8TinyVector32(
            data, work, first_inverse_skew, second_inverse_skew,
            forward_skew, 0);
        if (byte_count > 32)
            AVX2FF8HighEncodeTwoBlocksT8TinyVector32(
                data, work, first_inverse_skew, second_inverse_skew,
                forward_skew, byte_count - 32);
        return;
    }

    alignas(32) uint8_t staged[16][32] = {};
    const void* staged_data[16];
    void* staged_work[8];
    for (unsigned lane = 0; lane < 16; ++lane)
    {
        std::memcpy(
            staged[lane],
            static_cast<const uint8_t*>(data[lane]), byte_count);
        staged_data[lane] = staged[lane];
        if (lane < 8)
            staged_work[lane] = staged[lane];
    }

    /*
        The vector helper consumes every first-block input before overwriting
        staged rows 0..7 with inverse coefficients.  Rows 8..15 remain intact
        for the second inverse transform.  The fused coefficient XOR and one
        forward transform then leave exact parity rows in staged[0..7].
    */
    AVX2FF8HighEncodeTwoBlocksT8TinyVector32(
        staged_data, staged_work,
        first_inverse_skew, second_inverse_skew, forward_skew, 0);
    for (unsigned lane = 0; lane < 8; ++lane)
        std::memcpy(
            static_cast<uint8_t*>(work[lane]), staged[lane], byte_count);
}

bool TryAVX2FF8HighEncodeT8TailB256Packed(
    const void* const* data,
    void* const* recovery,
    uint32_t original_count,
    uint32_t recovery_count)
{
    if (!data || !recovery || original_count < 10 || original_count > 11 ||
        recovery_count < 5 || recovery_count > 8)
        return false;
    if (original_count == 10 && recovery_count == 8)
    {
        AVX2FF8HighEncodeK10R8T8B256Pairwise(data, recovery);
        return true;
    }
#define LEO2_T8_TAIL_DISPATCH(TailCount) \
    switch (recovery_count) \
    { \
    case 5: \
        AVX2FF8HighEncodeT8TailB256<TailCount, 5>(data, recovery); \
        break; \
    case 6: \
        AVX2FF8HighEncodeT8TailB256<TailCount, 6>(data, recovery); \
        break; \
    case 7: \
        AVX2FF8HighEncodeT8TailB256<TailCount, 7>(data, recovery); \
        break; \
    case 8: \
        AVX2FF8HighEncodeT8TailB256<TailCount, 8>(data, recovery); \
        break; \
    default: \
        return false; \
    }
    switch (original_count)
    {
    case 10:
        LEO2_T8_TAIL_DISPATCH(2);
        break;
    case 11:
        LEO2_T8_TAIL_DISPATCH(3);
        break;
    default:
        return false;
    }
#undef LEO2_T8_TAIL_DISPATCH
    return true;
}
#endif

#undef LEO2_AVX2_T8_ENTRY

#endif // promoted T=8 callback in the regular AVX2 backend

#if !defined(LEO2_AVX512_VARIANT)
/*
    The K=1/R=1 terminal has already proved that these ranges are disjoint.
    Keep eight cache lines in flight so the measured 4 KiB copy does not
    pay the libc dispatch/size-selection path.  Arbitrary tails remain correct
    for backend self-tests and future measured users, although production K1
    currently selects this callback only for 4096 bytes.
*/
static void AVX2CopyMemory(
    void* destination,
    const void* source,
    uint64_t byte_count)
{
    uint8_t* output = static_cast<uint8_t*>(destination);
    const uint8_t* input = static_cast<const uint8_t*>(source);
    while (byte_count >= 256)
    {
        const __m256i value0 = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(input));
        const __m256i value1 = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(input + 32));
        const __m256i value2 = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(input + 64));
        const __m256i value3 = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(input + 96));
        const __m256i value4 = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(input + 128));
        const __m256i value5 = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(input + 160));
        const __m256i value6 = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(input + 192));
        const __m256i value7 = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(input + 224));
        _mm256_storeu_si256(reinterpret_cast<__m256i*>(output), value0);
        _mm256_storeu_si256(
            reinterpret_cast<__m256i*>(output + 32), value1);
        _mm256_storeu_si256(
            reinterpret_cast<__m256i*>(output + 64), value2);
        _mm256_storeu_si256(
            reinterpret_cast<__m256i*>(output + 96), value3);
        _mm256_storeu_si256(
            reinterpret_cast<__m256i*>(output + 128), value4);
        _mm256_storeu_si256(
            reinterpret_cast<__m256i*>(output + 160), value5);
        _mm256_storeu_si256(
            reinterpret_cast<__m256i*>(output + 192), value6);
        _mm256_storeu_si256(
            reinterpret_cast<__m256i*>(output + 224), value7);
        input += 256;
        output += 256;
        byte_count -= 256;
    }
    while (byte_count >= 32)
    {
        const __m256i value = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(input));
        _mm256_storeu_si256(reinterpret_cast<__m256i*>(output), value);
        input += 32;
        output += 32;
        byte_count -= 32;
    }
    if (byte_count != 0)
        std::memcpy(output, input, static_cast<size_t>(byte_count));
}
#endif

static const Ops AVX2Ops = {
    LEO2_AVX_BACKEND_KIND,
    LEO2_AVX_BACKEND_NAME,
#ifdef LEO_HAS_FF8
    AVX2FF8Multiply,
    AVX2FF8MultiplyAdd,
#else
    NULL,
    NULL,
#endif
#ifdef LEO_HAS_FF16
    AVX2FF16Multiply,
    AVX2FF16MultiplyAdd,
#else
    NULL,
    NULL,
#endif
    AVX2XorMemory,
    AVX2XorMemory2To1,
    AVX2XorMemorySources,
    AVX2XorMemory4,
#ifdef LEO_HAS_FF8
    AVX2FF8IFFTButterfly2,
    AVX2FF8FFTButterfly2,
    AVX2FF8FFTButterfly2Out,
    AVX2FF8IFFTButterfly2Xor,
    AVX2FF8IFFTButterfly4,
    AVX2FF8FFTButterfly4,
    AVX2FF8IFFTButterfly4Out,
    AVX2FF8WeightedIFFTButterfly4,
    AVX2FF8FFTButterfly4Out,
    AVX2FF8IFFTButterfly4Range,
    AVX2FF8FFTButterfly4Range,
    AVX2FF8IFFTButterfly4XorRange,
#else
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
#endif
#ifdef LEO_HAS_FF16
    AVX2FF16IFFTButterfly2,
    AVX2FF16FFTButterfly2,
    AVX2FF16FFTButterfly2Out,
    AVX2FF16IFFTButterfly2Xor,
    AVX2FF16IFFTButterfly4,
    AVX2FF16FFTButterfly4,
    AVX2FF16IFFTButterfly4Out,
    AVX2FF16FFTButterfly4Out,
    AVX2FF16IFFTButterfly4Range,
    AVX2FF16FFTButterfly4Range
#else
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL
#endif
#if defined(LEO_HAS_FF8) && defined(LEO2_AVX512_VARIANT)
    , AVX2FF8HighEncodeOneBlock
#elif defined(LEO_HAS_FF8) && \
      !defined(LEO2_DIAGNOSTIC_DISABLE_HIGH_T8_VECTOR) && \
      !defined(LEO2_GFNI_VARIANT)
    , AVX2FF8HighEncodeOneBlockT8Vector
#else
    , NULL
#endif
#if defined(LEO_HAS_FF8) && defined(LEO2_AVX512_VARIANT)
    , kFF8HighEncodeSupportedSides
#elif defined(LEO_HAS_FF8) && \
      !defined(LEO2_DIAGNOSTIC_DISABLE_HIGH_T8_VECTOR) && \
      !defined(LEO2_GFNI_VARIANT)
    , 8U
#else
    , 0
#endif
#if defined(LEO_HAS_FF8) && \
    LEO2_EXPERIMENT_HIGH_T8_TWO_BLOCK_BINDING && \
    !defined(LEO2_DIAGNOSTIC_DISABLE_HIGH_T8_VECTOR) && \
    !defined(LEO2_AVX512_VARIANT) && !defined(LEO2_GFNI_VARIANT)
    , AVX2FF8HighEncodeTwoBlocksT8
#else
    , NULL
#endif
#if defined(LEO_HAS_FF8) && \
    LEO2_EXPERIMENT_HIGH_T8_TWO_BLOCK_BINDING && \
    !defined(LEO2_DIAGNOSTIC_DISABLE_HIGH_T8_VECTOR) && \
    !defined(LEO2_AVX512_VARIANT) && !defined(LEO2_GFNI_VARIANT)
    , AVX2FF8HighEncodeTwoBlocksT8Tiny
#else
    , NULL
#endif
#if defined(LEO_HAS_FF8) && !defined(LEO2_AVX512_VARIANT)
    , AVX2FF8HighEncodeSmall
#else
    , NULL
#endif
#if !defined(LEO2_AVX512_VARIANT)
    , AVX2XorMemoryDense
#else
    , NULL
#endif
#if defined(LEO_HAS_FF8) && !defined(LEO2_AVX512_VARIANT)
    , AVX2FF8MultiplyAddOutputs
#else
    , NULL
#endif
#if defined(LEO_HAS_FF8) && !defined(LEO2_AVX512_VARIANT)
    , AVX2FF8MultiplyAdd2Sources2Outputs
#else
    , NULL
#endif
#if defined(LEO_HAS_FF8) && LEO2_EXPERIMENT_GENERAL_ONE_LOSS_DIRECT && \
    !defined(LEO2_AVX512_VARIANT) && !defined(LEO2_GFNI_VARIANT)
    , AVX2FF8LinearCombination2
#else
    , NULL
#endif
#if !defined(LEO2_AVX512_VARIANT) && defined(LEO_HAS_FF8)
    , AVX2XorMemorySourcesFusedFinal
#else
    , NULL
#endif
#if defined(LEO2_GFNI_VARIANT) && defined(LEO_HAS_FF8)
    , AVX2FF8IFFTButterfly8Out
    , AVX2FF8FFTButterfly8Out
#else
    , NULL
    , NULL
#endif
#if defined(LEO_HAS_FF8) && LEO2_EXPERIMENT_GENERAL_ONE_LOSS_DIRECT && \
    !defined(LEO2_AVX512_VARIANT) && !defined(LEO2_GFNI_VARIANT)
    , AVX2FF8LinearCombination4Tiny
#else
    , NULL
#endif
#if defined(LEO_HAS_FF8) && !defined(LEO2_AVX512_VARIANT) && \
    !defined(LEO2_GFNI_VARIANT)
    , AVX2FF8HighEncodeT4Batch
#else
    , NULL
#endif
#if !defined(LEO2_AVX512_VARIANT)
    , AVX2CopyMemory
#else
    , NULL
#endif
#if defined(LEO_HAS_FF8) && !defined(LEO2_AVX512_VARIANT) && \
    !defined(LEO2_GFNI_VARIANT)
    , AVX2XorMemorySourcesGroup4
#else
    , NULL
#endif
#if defined(LEO_HAS_FF8) && !defined(LEO2_AVX512_VARIANT) && \
    !defined(LEO2_GFNI_VARIANT)
    , AVX2FF8IFFTButterfly2Range
#else
    , NULL
#endif
#if defined(LEO_HAS_FF8) && !defined(LEO2_AVX512_VARIANT) && \
    !defined(LEO2_GFNI_VARIANT)
    , AVX2FF8WalshLocator
#else
    , NULL
#endif
#if defined(LEO_HAS_FF8) && !defined(LEO2_AVX512_VARIANT) && \
    !defined(LEO2_GFNI_VARIANT)
    , AVX2XorMemorySourcesFixed64
    , AVX2XorMemorySourcesFixed256
#else
    , NULL
    , NULL
#endif
};

#if defined(LEO2_GFNI_VARIANT)
// Derives the affine operands from the same multiplication callbacks that build
// the nibble tables, so both backends are generated from one field definition.
// Fixed multiplication is GF(2)-linear, so each output byte of the product is
// determined by the products of the input basis vectors.
static bool GFNITablesPublished = false;

static const Ops* GFNIPublishTables(const InitializeArgs& args)
{
    if (GFNITablesPublished)
        return &AVX2Ops;
# ifdef LEO_HAS_FF8
#  ifdef LEO2_ENABLE_TEST_HOOKS
    if (TestShouldFailAllocation(LEO2_AVX_BACKEND_KIND, false))
        return NULL;
#  endif
    std::unique_ptr<FF8NibbleTable[]> ff8(
        new (std::nothrow) FF8NibbleTable[256]);
    if (!ff8)
        return NULL;
    for (unsigned log = 0; log < 256; ++log)
    {
        uint64_t matrix = 0;
        for (unsigned input_bit = 0; input_bit < 8; ++input_bit)
        {
            const uint8_t product = args.ff8_multiply_log(
                static_cast<uint8_t>(1U << input_bit),
                static_cast<uint8_t>(log));
            for (unsigned output_bit = 0; output_bit < 8; ++output_bit)
                if ((product >> output_bit) & 1U)
                    matrix |= GFNIAffineMatrixBit(output_bit, input_bit);
        }
        GFNIStoreMatrix(ff8[log].low, matrix);
        GFNIStoreMatrix(ff8[log].high, matrix);
    }
# endif
# ifdef LEO_HAS_FF16
#  ifdef LEO2_ENABLE_TEST_HOOKS
    if (TestShouldFailAllocation(LEO2_AVX_BACKEND_KIND, true))
        return NULL;
#  endif
    std::unique_ptr<FF16AffineTable[]> ff16(
        new (std::nothrow) FF16AffineTable[65536]);
    if (!ff16)
        return NULL;
    for (unsigned log = 0; log < 65536; ++log)
    {
        // block order: (low<-low, low<-high, high<-low, high<-high)
        uint64_t block[4] = { 0, 0, 0, 0 };
        for (unsigned input_bit = 0; input_bit < 8; ++input_bit)
        {
            const uint16_t from_low = args.ff16_multiply_log(
                static_cast<uint16_t>(1U << input_bit),
                static_cast<uint16_t>(log));
            const uint16_t from_high = args.ff16_multiply_log(
                static_cast<uint16_t>(1U << (input_bit + 8)),
                static_cast<uint16_t>(log));
            for (unsigned output_bit = 0; output_bit < 8; ++output_bit)
            {
                const uint64_t bit =
                    GFNIAffineMatrixBit(output_bit, input_bit);
                if ((from_low >> output_bit) & 1U)
                    block[0] |= bit;
                if ((from_high >> output_bit) & 1U)
                    block[1] |= bit;
                if ((from_low >> (output_bit + 8)) & 1U)
                    block[2] |= bit;
                if ((from_high >> (output_bit + 8)) & 1U)
                    block[3] |= bit;
            }
        }
        for (unsigned index = 0; index < 4; ++index)
            ff16[log].block[index] = block[index];
    }
# endif
# ifdef LEO_HAS_FF8
    FF8Tables = ff8.release();
# endif
# ifdef LEO_HAS_FF16
    FF16Tables = ff16.release();
# endif
    GFNITablesPublished = true;
    return &AVX2Ops;
}
#endif // LEO2_GFNI_VARIANT

const Ops* LEO2_AVX_INITIALIZER(const InitializeArgs& args)
{
#ifdef LEO_HAS_FF8
    if (!args.ff8_multiply_log)
        return NULL;
#endif
#ifdef LEO_HAS_FF16
    if (!args.ff16_multiply_log)
        return NULL;
#endif
#if defined(LEO2_GFNI_VARIANT) && !defined(LEO2_AVX512_VARIANT)
    // GFNI compiled under the AVX2 backend identity: this translation unit owns
    // the published table and no nibble table is required.
    return GFNIPublishTables(args);
#endif
#if defined(LEO2_AVX512_VARIANT)
# ifdef LEO2_ENABLE_TEST_HOOKS
#  ifdef LEO_HAS_FF8
    if (TestShouldFailAllocation(LEO2_AVX_BACKEND_KIND, false))
        return NULL;
#  endif
#  ifdef LEO_HAS_FF16
    if (TestShouldFailAllocation(LEO2_AVX_BACKEND_KIND, true))
        return NULL;
#  endif
# endif
    const Ops* const base = InitializeAVX2(args);
    if (!base)
        return NULL;
# if defined(LEO2_GFNI_VARIANT)
    return GFNIPublishTables(args);
# else
# ifdef LEO_HAS_FF8
    FF8Tables = static_cast<const FF8NibbleTable*>(GetAVX2FF8Tables());
    if (!FF8Tables)
        return NULL;
# endif
# ifdef LEO_HAS_FF16
    FF16Tables = static_cast<const FF16NibbleTable*>(GetAVX2FF16Tables());
    if (!FF16Tables)
        return NULL;
# endif
    return &AVX2Ops;
# endif
#else
#if defined(LEO_HAS_FF8) && defined(LEO_HAS_FF16)
    if (FF8Tables || FF16Tables)
        return FF8Tables && FF16Tables ? &AVX2Ops : NULL;
#elif defined(LEO_HAS_FF8)
    if (FF8Tables)
        return &AVX2Ops;
#else
    if (FF16Tables)
        return &AVX2Ops;
#endif
#ifdef LEO_HAS_FF8
#ifdef LEO2_ENABLE_TEST_HOOKS
    if (TestShouldFailAllocation(LEO2_AVX_BACKEND_KIND, false))
        return NULL;
#endif
    std::unique_ptr<FF8NibbleTable[]> ff8(
        new (std::nothrow) FF8NibbleTable[256]);
    if (!ff8)
        return NULL;
#endif
#if defined(LEO_HAS_FF16) && !defined(LEO2_GFNI_VARIANT)
#ifdef LEO2_ENABLE_TEST_HOOKS
    if (TestShouldFailAllocation(LEO2_AVX_BACKEND_KIND, true))
        return NULL;
#endif
    std::unique_ptr<FF16NibbleTable[]> ff16(
        new (std::nothrow) FF16NibbleTable[65536]);
    if (!ff16)
        return NULL;
#endif

#ifdef LEO_HAS_FF8
    for (unsigned log = 0; log < 256; ++log)
    {
        for (unsigned value = 0; value < 16; ++value)
        {
            ff8[log].low[value] = args.ff8_multiply_log(
                static_cast<uint8_t>(value), static_cast<uint8_t>(log));
            ff8[log].high[value] = args.ff8_multiply_log(
                static_cast<uint8_t>(value << 4), static_cast<uint8_t>(log));
        }
    }
#endif

#if defined(LEO_HAS_FF16) && !defined(LEO2_GFNI_VARIANT)
    for (int log = 0; log < 65536; ++log)
    {
        for (unsigned nibble = 0; nibble < 4; ++nibble)
        {
            for (unsigned value = 0; value < 16; ++value)
            {
                const uint16_t product = args.ff16_multiply_log(
                    static_cast<uint16_t>(value << (nibble * 4)),
                    static_cast<uint16_t>(log));
                ff16[log].low[nibble][value] =
                    static_cast<uint8_t>(product);
                ff16[log].high[nibble][value] =
                    static_cast<uint8_t>(product >> 8);
            }
        }
    }
#endif
#ifdef LEO_HAS_FF8
    FF8Tables = ff8.release();
#endif
#if defined(LEO_HAS_FF16) && !defined(LEO2_GFNI_VARIANT)
    FF16Tables = ff16.release();
#endif
    return &AVX2Ops;
#endif
}

#if !defined(LEO2_AVX512_VARIANT) && !defined(LEO2_GFNI_MEMBER)
const void* GetAVX2FF8Tables()
{
#ifdef LEO_HAS_FF8
    return FF8Tables;
#else
    return NULL;
#endif
}

const void* GetAVX2FF16Tables()
{
#ifdef LEO_HAS_FF16
    return FF16Tables;
#else
    return NULL;
#endif
}
#endif

/*
    The ordinary AVX2 member and the documented in-place GFNI diagnostic both
    own these concrete symbols.  The separately linked production GFNI member
    must not emit a second copy; its core dispatch never calls them under the
    GFNI backend identity.
*/
#if defined(LEO_HAS_FF8) && !defined(LEO2_AVX512_VARIANT) && \
    !defined(LEO2_GFNI_MEMBER)
template<uint32_t OriginalCount>
static LEO_FORCE_INLINE void AVX2FF8HighEncodeT2MultiVector(
    const uint8_t* const* input,
    uint8_t* output0,
    uint8_t* output1,
    uint64_t offset,
    const AVX2FF8LinearTable* multipliers)
{
    static_assert(OriginalCount >= 5 && OriginalCount <= 16,
        "T=2 multi-block circuit instantiated outside its K range");
    __m256i parity0 = _mm256_setzero_si256();
    __m256i parity1 = _mm256_setzero_si256();
    uint32_t source = 0;
    uint32_t pair = 0;
    for (; source + 1U < OriginalCount; source += 2U, ++pair)
    {
        const __m256i a = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(input[source] + offset));
        const __m256i b = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(input[source + 1U] + offset));
        const __m256i common = AVX2FF8ProductVector(
            _mm256_xor_si256(a, b),
            multipliers[pair].low, multipliers[pair].high);
        parity0 = _mm256_xor_si256(
            parity0, _mm256_xor_si256(a, common));
        parity1 = _mm256_xor_si256(
            parity1, _mm256_xor_si256(b, common));
    }
    if (source < OriginalCount)
    {
        const __m256i a = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(input[source] + offset));
        const __m256i common = AVX2FF8ProductVector(
            a, multipliers[pair].low, multipliers[pair].high);
        parity0 = _mm256_xor_si256(
            parity0, _mm256_xor_si256(a, common));
        parity1 = _mm256_xor_si256(parity1, common);
    }
    _mm256_storeu_si256(
        reinterpret_cast<__m256i*>(output0 + offset), parity0);
    _mm256_storeu_si256(
        reinterpret_cast<__m256i*>(output1 + offset), parity1);
}

template<uint32_t OriginalCount>
static void AVX2FF8HighEncodeT2MultiPackedShape(
    const void* const* data,
    void* const* recovery,
    uint64_t byte_count)
{
    static_assert(OriginalCount >= 5 && OriginalCount <= 16,
        "T=2 multi-block circuit instantiated outside its K range");
    static const uint16_t kMultiplierLogs[8] = {
        85, 17, 34, 153, 102, 51, 187, 219
    };
    static const uint16_t kZeroSkew = 255;
    const uint32_t pair_count = (OriginalCount + 1U) / 2U;
    AVX2FF8LinearTable multipliers[8];
    for (uint32_t pair = 0; pair < pair_count; ++pair)
    {
        multipliers[pair] = AVX2FF8PrepareLinearTable(
            kMultiplierLogs[pair], kZeroSkew);
    }
    const uint8_t* input[OriginalCount];
    for (uint32_t i = 0; i < OriginalCount; ++i)
        input[i] = static_cast<const uint8_t*>(data[i]);
    uint8_t* const output0 = static_cast<uint8_t*>(recovery[0]);
    uint8_t* const output1 = static_cast<uint8_t*>(recovery[1]);
    for (uint64_t offset = 0; offset < byte_count; offset += 32U)
    {
        AVX2FF8HighEncodeT2MultiVector<OriginalCount>(
            input, output0, output1, offset, multipliers);
    }
}

void AVX2FF8HighEncodeT2MultiPacked(
    const void* const* data,
    void* const* recovery,
    uint32_t original_count,
    uint64_t byte_count)
{
    LEO_DEBUG_ASSERT(data && recovery && original_count >= 5 &&
        original_count <= 16 && byte_count >= 64U &&
        (byte_count & 63U) == 0);
#define LEO2_T2_MULTI_CASE(K) \
    case K: AVX2FF8HighEncodeT2MultiPackedShape<K>( \
        data, recovery, byte_count); return
    switch (original_count)
    {
    LEO2_T2_MULTI_CASE(5);
    LEO2_T2_MULTI_CASE(6);
    LEO2_T2_MULTI_CASE(7);
    LEO2_T2_MULTI_CASE(8);
    LEO2_T2_MULTI_CASE(9);
    LEO2_T2_MULTI_CASE(10);
    LEO2_T2_MULTI_CASE(11);
    LEO2_T2_MULTI_CASE(12);
    LEO2_T2_MULTI_CASE(13);
    LEO2_T2_MULTI_CASE(14);
    LEO2_T2_MULTI_CASE(15);
    LEO2_T2_MULTI_CASE(16);
    default:
        LEO_DEBUG_BREAK;
        return;
    }
#undef LEO2_T2_MULTI_CASE
}

void AVX2FF8HighEncodeT2Packed(
    const void* const* data,
    void* const* recovery,
    uint32_t original_count,
    uint64_t byte_count)
{
    LEO_DEBUG_ASSERT(data != NULL && recovery != NULL);
    LEO_DEBUG_ASSERT(original_count == 2 || original_count == 3);
    LEO_DEBUG_ASSERT(byte_count >= 32 && (byte_count & 31U) == 0);

    /*
        In the legacy Cantor representation, log(2)=85 and log(4)=17.
        Expanding the fixed T=2 inverse/forward circuit gives

          p0 = 3*a + 2*b,     p1 = 2*a + 3*b                  (K=2)
          p0 ^= 5*c,          p1 ^= 4*c                       (K=3)

        or equivalently the existing common-term form with multipliers 2 and
        4.  The direct algebra tests independently verify these rows against
        the systematic generator matrix and the padded transform.
    */
    static const uint16_t kZeroSkew = 255;
    const AVX2FF8LinearTable first =
        AVX2FF8PrepareLinearTable(85, kZeroSkew);
    const uint8_t* input[3] = {
        static_cast<const uint8_t*>(data[0]),
        static_cast<const uint8_t*>(data[1]),
        original_count == 3
            ? static_cast<const uint8_t*>(data[2])
            : NULL
    };
    uint8_t* output0 = static_cast<uint8_t*>(recovery[0]);
    uint8_t* output1 = static_cast<uint8_t*>(recovery[1]);

    if (original_count == 2)
    {
        for (uint64_t offset = 0; offset < byte_count; offset += 32U)
        {
            AVX2FF8HighEncodeT2Vector<2>(
                input, output0, output1, offset, first, first);
        }
        return;
    }
    const AVX2FF8LinearTable second =
        AVX2FF8PrepareLinearTable(17, kZeroSkew);
    for (uint64_t offset = 0; offset < byte_count; offset += 32U)
    {
        AVX2FF8HighEncodeT2Vector<3>(
            input, output0, output1, offset, first, second);
    }
}

void AVX2FF8HighEncodeT2PackedTail(
    const void* const* data,
    void* const* recovery,
    uint32_t original_count,
    uint64_t byte_count)
{
    LEO_DEBUG_ASSERT(data != NULL && recovery != NULL);
    LEO_DEBUG_ASSERT(original_count == 2 || original_count == 3);
    LEO_DEBUG_ASSERT(byte_count != 0 && (byte_count & 31U) != 0);

    const uint8_t* input[3] = {
        static_cast<const uint8_t*>(data[0]),
        static_cast<const uint8_t*>(data[1]),
        original_count == 3
            ? static_cast<const uint8_t*>(data[2])
            : NULL
    };
    uint8_t* output0 = static_cast<uint8_t*>(recovery[0]);
    uint8_t* output1 = static_cast<uint8_t*>(recovery[1]);
    const uint64_t aligned_bytes = byte_count & ~UINT64_C(31);
    if (aligned_bytes != 0)
    {
        AVX2FF8HighEncodeT2Packed(
            data, recovery, original_count, aligned_bytes);
    }

    /* The table representation is the scalar oracle for the same two fixed
       maps used above.  Evaluating only the final ragged bytes avoids padded
       loads while retaining the one-pass circuit. */
    for (uint64_t offset = aligned_bytes; offset < byte_count; ++offset)
    {
        const uint8_t a = input[0][offset];
        const uint8_t b = input[1][offset];
        uint8_t common = FF8Product(85,
            static_cast<uint8_t>(a ^ b));
        uint8_t parity0 = a;
        if (original_count == 3)
        {
            const uint8_t c = input[2][offset];
            common = static_cast<uint8_t>(common ^ FF8Product(17, c));
            parity0 = static_cast<uint8_t>(parity0 ^ c);
        }
        output0[offset] = static_cast<uint8_t>(parity0 ^ common);
        output1[offset] = static_cast<uint8_t>(b ^ common);
    }
}
#endif

#ifdef LEO2_ENABLE_TEST_HOOKS
void LEO2_AVX_TABLE_STATE(TestBackendState* state)
{
#ifdef LEO_HAS_FF8
    state->ff8_published = FF8Tables != NULL;
    state->ff8_bytes = 256U * sizeof(FF8NibbleTable);
#else
    state->ff8_published = false;
    state->ff8_bytes = 0;
#endif
#ifdef LEO_HAS_FF16
    state->ff16_published = FF16Tables != NULL;
    state->ff16_bytes = 65536U * sizeof(FF16Table);
#else
    state->ff16_published = false;
    state->ff16_bytes = 0;
#endif
}
#endif

#undef LEO2_AVX_TABLE_STATE
#undef LEO2_AVX_INITIALIZER
#undef LEO2_AVX_BACKEND_NAME
#undef LEO2_AVX_BACKEND_KIND

}} // namespace leopard::backend
