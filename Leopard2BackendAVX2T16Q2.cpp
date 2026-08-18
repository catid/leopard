/*
    Copyright (c) 2017 Christopher A. Taylor.  All rights reserved.
    See Leopard2Backend.h for the full BSD license notice.
*/

#include "Leopard2Backend.h"
#include "LeopardCommon.h"

#include <immintrin.h>

namespace leopard { namespace backend {

#if LEO2_EXPERIMENT_HIGH_T16_Q2_B64_FUSED && !defined(NO_LEO_HAS_FF8)

/* Published by the ordinary AVX2 translation unit after table setup. */
const void* GetAVX2FF8Tables();

/*
    The legacy-high encoder advances the inverse-transform skew by T for each
    complete message block.  For T=16, block b therefore reads the 16 factors
    at FFTSkewStorage + (b + 1) * 16, while parity evaluation reads offset 0:

      forward +0:
        0,255,255,85,255,17,85,34,255,153,17,102,85,51,34,187
      inverse block 0, +16:
        255,219,153,7,17,111,102,28,85,183,51,224,34,131,187,222
      inverse block 1, +32:
        255,196,219,76,153,54,7,99,17,19,111,49,102,67,28,52

    Radix-four group r consumes [r+1], [r+3], [r+2].  Its outer distance-four
    group consumes [4], [12], [8].  These index rules produce every template
    argument below and match IFFT_DIT_Encoder/FFT_DIT exactly.
*/

static const size_t kT16Q2TableBytes = 32;

static LEO_FORCE_INLINE __m256i T16Q2Broadcast(const uint8_t table[16])
{
    return _mm256_broadcastsi128_si256(_mm_loadu_si128(
        reinterpret_cast<const __m128i*>(table)));
}

static LEO_FORCE_INLINE __m256i T16Q2Product(
    __m256i data,
    const unsigned char* table)
{
    const __m256i nibble_mask = _mm256_set1_epi8(15);
    const __m256i low = _mm256_shuffle_epi8(
        T16Q2Broadcast(table),
        _mm256_and_si256(data, nibble_mask));
    const __m256i high = _mm256_shuffle_epi8(
        T16Q2Broadcast(table + 16),
        _mm256_and_si256(_mm256_srli_epi64(data, 4), nibble_mask));
    __m256i product = _mm256_xor_si256(low, high);
#if defined(__GNUC__) || defined(__clang__)
    /*
        Keep each constant-table pair local to this product.  Without this
        compiler barrier GCC hoists several pairs across the composed
        eight-input transform, exhausting the sixteen AVX2 registers and
        spilling live codeword vectors to the stack.  The register operand
        also prevents GCC from splitting one product around the next stage.
    */
    __asm__ __volatile__("" : "+x"(product) : : "memory");
#endif
    return product;
}

static LEO_FORCE_INLINE __m256i T16Q2Load(
    const void* pointer,
    uint64_t offset)
{
    const uint8_t* bytes = static_cast<const uint8_t*>(pointer) + offset;
    return _mm256_loadu_si256(reinterpret_cast<const __m256i*>(bytes));
}

static LEO_FORCE_INLINE void T16Q2Store(
    void* pointer,
    uint64_t offset,
    __m256i value)
{
    uint8_t* bytes = static_cast<uint8_t*>(pointer) + offset;
    _mm256_storeu_si256(reinterpret_cast<__m256i*>(bytes), value);
}

template<unsigned Log>
static LEO_FORCE_INLINE void T16Q2MulAdd(
    const unsigned char* tables,
    __m256i& destination,
    __m256i source)
{
    static_assert(Log <= 255, "T16 fixed multiplier log is out of range");
    destination = _mm256_xor_si256(
        destination, T16Q2Product(
            source, tables + Log * kT16Q2TableBytes));
}

template<unsigned Log01, unsigned Log23, unsigned Log02>
static LEO_FORCE_INLINE void T16Q2IFFT4(
    const unsigned char* tables,
    __m256i& value0,
    __m256i& value1,
    __m256i& value2,
    __m256i& value3)
{
    value1 = _mm256_xor_si256(value1, value0);
    T16Q2MulAdd<Log01>(tables, value0, value1);
    value3 = _mm256_xor_si256(value3, value2);
    T16Q2MulAdd<Log23>(tables, value2, value3);
    value2 = _mm256_xor_si256(value2, value0);
    value3 = _mm256_xor_si256(value3, value1);
    T16Q2MulAdd<Log02>(tables, value0, value2);
    T16Q2MulAdd<Log02>(tables, value1, value3);
}

template<unsigned Log01, unsigned Log23, unsigned Log02>
static LEO_FORCE_INLINE void T16Q2FFT4(
    const unsigned char* tables,
    __m256i& value0,
    __m256i& value1,
    __m256i& value2,
    __m256i& value3)
{
    if (Log02 != 255)
    {
        T16Q2MulAdd<Log02>(tables, value0, value2);
        T16Q2MulAdd<Log02>(tables, value1, value3);
    }
    value2 = _mm256_xor_si256(value2, value0);
    value3 = _mm256_xor_si256(value3, value1);
    if (Log01 != 255)
        T16Q2MulAdd<Log01>(tables, value0, value1);
    value1 = _mm256_xor_si256(value1, value0);
    if (Log23 != 255)
        T16Q2MulAdd<Log23>(tables, value2, value3);
    value3 = _mm256_xor_si256(value3, value2);
}

template<
    uint64_t ShardBytes,
    unsigned Log01,
    unsigned Log23,
    unsigned Log02,
    unsigned InputRow,
    unsigned OutputRow,
    bool PartialInput>
static LEO_FORCE_INLINE void T16Q2InverseGroup(
    const unsigned char* tables,
    const void* data_base,
    void* work_base,
    unsigned original_count)
{
    if (PartialInput && InputRow >= original_count)
    {
        const __m256i zero = _mm256_setzero_si256();
        for (uint64_t offset = 0; offset < ShardBytes; offset += 32)
        {
            T16Q2Store(work_base,
                (OutputRow + 0U) * ShardBytes + offset, zero);
            T16Q2Store(work_base,
                (OutputRow + 1U) * ShardBytes + offset, zero);
            T16Q2Store(work_base,
                (OutputRow + 2U) * ShardBytes + offset, zero);
            T16Q2Store(work_base,
                (OutputRow + 3U) * ShardBytes + offset, zero);
        }
        return;
    }

    for (uint64_t offset = 0; offset < ShardBytes; offset += 32)
    {
        const __m256i zero = _mm256_setzero_si256();
        __m256i value0 = !PartialInput || InputRow + 0U < original_count
            ? T16Q2Load(data_base,
                (InputRow + 0U) * ShardBytes + offset)
            : zero;
        __m256i value1 = !PartialInput || InputRow + 1U < original_count
            ? T16Q2Load(data_base,
                (InputRow + 1U) * ShardBytes + offset)
            : zero;
        __m256i value2 = !PartialInput || InputRow + 2U < original_count
            ? T16Q2Load(data_base,
                (InputRow + 2U) * ShardBytes + offset)
            : zero;
        __m256i value3 = !PartialInput || InputRow + 3U < original_count
            ? T16Q2Load(data_base,
                (InputRow + 3U) * ShardBytes + offset)
            : zero;
        T16Q2IFFT4<Log01, Log23, Log02>(
            tables, value0, value1, value2, value3);
        T16Q2Store(work_base,
            (OutputRow + 0U) * ShardBytes + offset, value0);
        T16Q2Store(work_base,
            (OutputRow + 1U) * ShardBytes + offset, value1);
        T16Q2Store(work_base,
            (OutputRow + 2U) * ShardBytes + offset, value2);
        T16Q2Store(work_base,
            (OutputRow + 3U) * ShardBytes + offset, value3);
    }
}

template<uint64_t ShardBytes, unsigned ActiveSecondGroups>
static LEO_FORCE_INLINE void T16Q2FusedOuter(
    const unsigned char* tables,
    void* recovery_base,
    void* temporary_base)
{
    for (unsigned column = 0; column < 4; ++column)
    {
        for (uint64_t offset = 0; offset < ShardBytes; offset += 32)
        {
            const __m256i zero = _mm256_setzero_si256();
            __m256i first0 = T16Q2Load(
                recovery_base,
                (column + 0U) * ShardBytes + offset);
            __m256i first1 = T16Q2Load(
                recovery_base,
                (column + 4U) * ShardBytes + offset);
            __m256i first2 = T16Q2Load(
                recovery_base,
                (column + 8U) * ShardBytes + offset);
            __m256i first3 = T16Q2Load(
                recovery_base,
                (column + 12U) * ShardBytes + offset);
            __m256i second0 = T16Q2Load(
                temporary_base,
                (column + 0U) * ShardBytes + offset);
            __m256i second1 = ActiveSecondGroups >= 2
                ? T16Q2Load(temporary_base,
                    (column + 4U) * ShardBytes + offset)
                : zero;
            __m256i second2 = ActiveSecondGroups >= 3
                ? T16Q2Load(temporary_base,
                    (column + 8U) * ShardBytes + offset)
                : zero;
            __m256i second3 = ActiveSecondGroups >= 4
                ? T16Q2Load(temporary_base,
                    (column + 12U) * ShardBytes + offset)
                : zero;

            /* Final IFFT layer for the first and second message cosets. */
            T16Q2IFFT4<17, 34, 85>(
                tables, first0, first1, first2, first3);
            T16Q2IFFT4<153, 102, 17>(
                tables, second0, second1, second2, second3);
            first0 = _mm256_xor_si256(first0, second0);
            first1 = _mm256_xor_si256(first1, second1);
            first2 = _mm256_xor_si256(first2, second2);
            first3 = _mm256_xor_si256(first3, second3);

            /* First forward layer on the accumulated coefficients. */
            T16Q2FFT4<255, 85, 255>(
                tables, first0, first1, first2, first3);

            T16Q2Store(
                recovery_base,
                (column + 0U) * ShardBytes + offset, first0);
            T16Q2Store(
                recovery_base,
                (column + 4U) * ShardBytes + offset, first1);
            T16Q2Store(
                recovery_base,
                (column + 8U) * ShardBytes + offset, first2);
            T16Q2Store(
                recovery_base,
                (column + 12U) * ShardBytes + offset, first3);
        }
    }
}

template<
    uint64_t ShardBytes,
    unsigned Log01,
    unsigned Log23,
    unsigned Log02,
    unsigned BaseRow,
    bool PartialOutput>
static LEO_FORCE_INLINE void T16Q2ForwardGroup(
    const unsigned char* tables,
    const void* work_base,
    void* recovery_base,
    unsigned recovery_count)
{
    if (PartialOutput && BaseRow >= recovery_count)
        return;
    for (uint64_t offset = 0; offset < ShardBytes; offset += 32)
    {
        __m256i value0 = T16Q2Load(
            work_base,
            (BaseRow + 0U) * ShardBytes + offset);
        __m256i value1 = T16Q2Load(
            work_base,
            (BaseRow + 1U) * ShardBytes + offset);
        __m256i value2 = T16Q2Load(
            work_base,
            (BaseRow + 2U) * ShardBytes + offset);
        __m256i value3 = T16Q2Load(
            work_base,
            (BaseRow + 3U) * ShardBytes + offset);
        T16Q2FFT4<Log01, Log23, Log02>(
            tables, value0, value1, value2, value3);
        if (!PartialOutput || BaseRow + 0U < recovery_count)
            T16Q2Store(recovery_base,
                (BaseRow + 0U) * ShardBytes + offset, value0);
        if (!PartialOutput || BaseRow + 1U < recovery_count)
            T16Q2Store(recovery_base,
                (BaseRow + 1U) * ShardBytes + offset, value1);
        if (!PartialOutput || BaseRow + 2U < recovery_count)
            T16Q2Store(recovery_base,
                (BaseRow + 2U) * ShardBytes + offset, value2);
        if (!PartialOutput || BaseRow + 3U < recovery_count)
            T16Q2Store(recovery_base,
                (BaseRow + 3U) * ShardBytes + offset, value3);
    }
}

template<
    uint64_t ShardBytes,
    unsigned ActiveSecondGroups,
    bool PartialInput,
    bool PartialOutput>
static LEO_FORCE_INLINE void T16Q2Encode(
    const unsigned char* tables,
    const void* data_base,
    void* recovery_base,
    void* temporary_base,
    unsigned original_count,
    unsigned recovery_count)
{
    uint8_t* const temporary_bytes =
        static_cast<uint8_t*>(temporary_base);
    void* const primary_base = PartialOutput
        ? temporary_base : recovery_base;
    void* const secondary_base = PartialOutput
        ? temporary_bytes + 16U * ShardBytes : temporary_base;

    /* First two inverse layers, FFTSkewStorage + 16. */
    T16Q2InverseGroup<ShardBytes, 219, 7, 153, 0, 0, false>(
        tables, data_base, primary_base, original_count);
    T16Q2InverseGroup<ShardBytes, 111, 28, 102, 4, 4, false>(
        tables, data_base, primary_base, original_count);
    T16Q2InverseGroup<ShardBytes, 183, 224, 51, 8, 8, false>(
        tables, data_base, primary_base, original_count);
    T16Q2InverseGroup<ShardBytes, 131, 222, 187, 12, 12, false>(
        tables, data_base, primary_base, original_count);

    /* First two inverse layers, FFTSkewStorage + 32. */
    T16Q2InverseGroup<ShardBytes, 196, 76, 219, 16, 0, PartialInput>(
        tables, data_base, secondary_base, original_count);
    if (ActiveSecondGroups >= 2)
        T16Q2InverseGroup<ShardBytes, 54, 99, 7, 20, 4, PartialInput>(
            tables, data_base, secondary_base, original_count);
    if (ActiveSecondGroups >= 3)
        T16Q2InverseGroup<ShardBytes, 19, 49, 111, 24, 8, PartialInput>(
            tables, data_base, secondary_base, original_count);
    if (ActiveSecondGroups >= 4)
        T16Q2InverseGroup<ShardBytes, 67, 52, 28, 28, 12, PartialInput>(
            tables, data_base, secondary_base, original_count);

    T16Q2FusedOuter<ShardBytes, ActiveSecondGroups>(
        tables, primary_base, secondary_base);

    /* Final forward layer, FFTSkewStorage + 0. */
    T16Q2ForwardGroup<ShardBytes, 255, 85, 255, 0, PartialOutput>(
        tables, primary_base, recovery_base, recovery_count);
    T16Q2ForwardGroup<ShardBytes, 17, 34, 85, 4, PartialOutput>(
        tables, primary_base, recovery_base, recovery_count);
    T16Q2ForwardGroup<ShardBytes, 153, 102, 17, 8, PartialOutput>(
        tables, primary_base, recovery_base, recovery_count);
    T16Q2ForwardGroup<ShardBytes, 51, 187, 34, 12, PartialOutput>(
        tables, primary_base, recovery_base, recovery_count);
}

template<uint64_t ShardBytes, unsigned ActiveSecondGroups>
static LEO_FORCE_INLINE void T16Q2EncodeSelected(
    const unsigned char* tables,
    const void* data_base,
    void* recovery_base,
    void* temporary_base,
    unsigned original_count,
    unsigned recovery_count)
{
    if (original_count == 32 && recovery_count == 16)
        T16Q2Encode<ShardBytes, ActiveSecondGroups, false, false>(
            tables, data_base,
            recovery_base, temporary_base, original_count, recovery_count);
    else if (recovery_count == 16)
        T16Q2Encode<ShardBytes, ActiveSecondGroups, true, false>(
            tables, data_base,
            recovery_base, temporary_base, original_count, recovery_count);
    else
        T16Q2Encode<ShardBytes, ActiveSecondGroups, true, true>(
            tables, data_base,
            recovery_base, temporary_base, original_count, recovery_count);
}

template<uint64_t ShardBytes, unsigned ActiveSecondGroups>
static LEO_FORCE_INLINE void T16Q2EncodePartialSelected(
    const unsigned char* tables,
    const void* data_base,
    void* recovery_base,
    void* temporary_base,
    unsigned original_count,
    unsigned recovery_count)
{
    if (recovery_count == 16)
        T16Q2Encode<ShardBytes, ActiveSecondGroups, true, false>(
            tables, data_base,
            recovery_base, temporary_base, original_count, recovery_count);
    else
        T16Q2Encode<ShardBytes, ActiveSecondGroups, true, true>(
            tables, data_base,
            recovery_base, temporary_base, original_count, recovery_count);
}

/*
    Extend the same transform to the third and fourth T=16 message cosets.
    The fixed tuples below are FFTSkewStorage slices +48 and +64, generated by
    Leopard's FFTInitialize() and independently checked against its +0/+16/+32
    slices:

      +48: 85,137,183,226,51,26,224,108,
           34,27,131,134,187,38,222,139
      +64: 255,241,196,254,219,31,76,239,
           153,73,54,200,7,148,99,140

    Keep the first pair's outer-IFFT result as a coefficient accumulator.
    Each later block reuses the second 16-row half, XORs its completed inverse
    transform into that accumulator, and the last block fuses the first
    forward layer.  This is the exact legacy block-IFFT sum followed by one
    parity FFT, with no additional shard workspace.
*/
static LEO_FORCE_INLINE void T16Q4MergeFirstPairOuter(
    const unsigned char* tables,
    void* primary_base,
    const void* secondary_base)
{
    for (unsigned column = 0; column < 4; ++column)
    {
        for (uint64_t offset = 0; offset < 64; offset += 32)
        {
            __m256i first0 = T16Q2Load(
                primary_base, (column + 0U) * 64U + offset);
            __m256i first1 = T16Q2Load(
                primary_base, (column + 4U) * 64U + offset);
            __m256i first2 = T16Q2Load(
                primary_base, (column + 8U) * 64U + offset);
            __m256i first3 = T16Q2Load(
                primary_base, (column + 12U) * 64U + offset);
            __m256i second0 = T16Q2Load(
                secondary_base, (column + 0U) * 64U + offset);
            __m256i second1 = T16Q2Load(
                secondary_base, (column + 4U) * 64U + offset);
            __m256i second2 = T16Q2Load(
                secondary_base, (column + 8U) * 64U + offset);
            __m256i second3 = T16Q2Load(
                secondary_base, (column + 12U) * 64U + offset);

            T16Q2IFFT4<17, 34, 85>(
                tables, first0, first1, first2, first3);
            T16Q2IFFT4<153, 102, 17>(
                tables, second0, second1, second2, second3);
            first0 = _mm256_xor_si256(first0, second0);
            first1 = _mm256_xor_si256(first1, second1);
            first2 = _mm256_xor_si256(first2, second2);
            first3 = _mm256_xor_si256(first3, second3);

            T16Q2Store(primary_base,
                (column + 0U) * 64U + offset, first0);
            T16Q2Store(primary_base,
                (column + 4U) * 64U + offset, first1);
            T16Q2Store(primary_base,
                (column + 8U) * 64U + offset, first2);
            T16Q2Store(primary_base,
                (column + 12U) * 64U + offset, first3);
        }
    }
}

template<unsigned Log01, unsigned Log23, unsigned Log02, bool FinishForward>
static LEO_FORCE_INLINE void T16Q4AccumulateOuter(
    const unsigned char* tables,
    void* primary_base,
    const void* secondary_base)
{
    for (unsigned column = 0; column < 4; ++column)
    {
        for (uint64_t offset = 0; offset < 64; offset += 32)
        {
            __m256i first0 = T16Q2Load(
                primary_base, (column + 0U) * 64U + offset);
            __m256i first1 = T16Q2Load(
                primary_base, (column + 4U) * 64U + offset);
            __m256i first2 = T16Q2Load(
                primary_base, (column + 8U) * 64U + offset);
            __m256i first3 = T16Q2Load(
                primary_base, (column + 12U) * 64U + offset);
            __m256i second0 = T16Q2Load(
                secondary_base, (column + 0U) * 64U + offset);
            __m256i second1 = T16Q2Load(
                secondary_base, (column + 4U) * 64U + offset);
            __m256i second2 = T16Q2Load(
                secondary_base, (column + 8U) * 64U + offset);
            __m256i second3 = T16Q2Load(
                secondary_base, (column + 12U) * 64U + offset);

            T16Q2IFFT4<Log01, Log23, Log02>(
                tables, second0, second1, second2, second3);
            first0 = _mm256_xor_si256(first0, second0);
            first1 = _mm256_xor_si256(first1, second1);
            first2 = _mm256_xor_si256(first2, second2);
            first3 = _mm256_xor_si256(first3, second3);
            if (FinishForward)
            {
                T16Q2FFT4<255, 85, 255>(
                    tables, first0, first1, first2, first3);
            }

            T16Q2Store(primary_base,
                (column + 0U) * 64U + offset, first0);
            T16Q2Store(primary_base,
                (column + 4U) * 64U + offset, first1);
            T16Q2Store(primary_base,
                (column + 8U) * 64U + offset, first2);
            T16Q2Store(primary_base,
                (column + 12U) * 64U + offset, first3);
        }
    }
}

static LEO_FORCE_INLINE void T16Q4Encode(
    const unsigned char* tables,
    const void* data_base,
    void* recovery_base,
    void* temporary_base,
    unsigned original_count,
    unsigned recovery_count)
{
    uint8_t* const temporary_bytes =
        static_cast<uint8_t*>(temporary_base);
    void* const primary_base = recovery_count == 16
        ? recovery_base : temporary_base;
    void* const secondary_base = recovery_count == 16
        ? temporary_base : temporary_bytes + 16U * 64U;

    T16Q2InverseGroup<64, 219, 7, 153, 0, 0, false>(
        tables, data_base, primary_base, original_count);
    T16Q2InverseGroup<64, 111, 28, 102, 4, 4, false>(
        tables, data_base, primary_base, original_count);
    T16Q2InverseGroup<64, 183, 224, 51, 8, 8, false>(
        tables, data_base, primary_base, original_count);
    T16Q2InverseGroup<64, 131, 222, 187, 12, 12, false>(
        tables, data_base, primary_base, original_count);
    T16Q2InverseGroup<64, 196, 76, 219, 16, 0, false>(
        tables, data_base, secondary_base, original_count);
    T16Q2InverseGroup<64, 54, 99, 7, 20, 4, false>(
        tables, data_base, secondary_base, original_count);
    T16Q2InverseGroup<64, 19, 49, 111, 24, 8, false>(
        tables, data_base, secondary_base, original_count);
    T16Q2InverseGroup<64, 67, 52, 28, 28, 12, false>(
        tables, data_base, secondary_base, original_count);
    T16Q4MergeFirstPairOuter(tables, primary_base, secondary_base);

    T16Q2InverseGroup<64, 137, 226, 183, 32, 0, true>(
        tables, data_base, secondary_base, original_count);
    T16Q2InverseGroup<64, 26, 108, 224, 36, 4, true>(
        tables, data_base, secondary_base, original_count);
    T16Q2InverseGroup<64, 27, 134, 131, 40, 8, true>(
        tables, data_base, secondary_base, original_count);
    T16Q2InverseGroup<64, 38, 139, 222, 44, 12, true>(
        tables, data_base, secondary_base, original_count);

    if (original_count <= 48)
    {
        T16Q4AccumulateOuter<51, 187, 34, true>(
            tables, primary_base, secondary_base);
    }
    else
    {
        T16Q4AccumulateOuter<51, 187, 34, false>(
            tables, primary_base, secondary_base);
        T16Q2InverseGroup<64, 241, 254, 196, 48, 0, true>(
            tables, data_base, secondary_base, original_count);
        T16Q2InverseGroup<64, 31, 239, 76, 52, 4, true>(
            tables, data_base, secondary_base, original_count);
        T16Q2InverseGroup<64, 73, 200, 54, 56, 8, true>(
            tables, data_base, secondary_base, original_count);
        T16Q2InverseGroup<64, 148, 140, 99, 60, 12, true>(
            tables, data_base, secondary_base, original_count);
        T16Q4AccumulateOuter<219, 7, 153, true>(
            tables, primary_base, secondary_base);
    }

    if (recovery_count == 16)
    {
        T16Q2ForwardGroup<64, 255, 85, 255, 0, false>(
            tables, primary_base, recovery_base, recovery_count);
        T16Q2ForwardGroup<64, 17, 34, 85, 4, false>(
            tables, primary_base, recovery_base, recovery_count);
        T16Q2ForwardGroup<64, 153, 102, 17, 8, false>(
            tables, primary_base, recovery_base, recovery_count);
        T16Q2ForwardGroup<64, 51, 187, 34, 12, false>(
            tables, primary_base, recovery_base, recovery_count);
    }
    else
    {
        T16Q2ForwardGroup<64, 255, 85, 255, 0, true>(
            tables, primary_base, recovery_base, recovery_count);
        T16Q2ForwardGroup<64, 17, 34, 85, 4, true>(
            tables, primary_base, recovery_base, recovery_count);
        T16Q2ForwardGroup<64, 153, 102, 17, 8, true>(
            tables, primary_base, recovery_base, recovery_count);
        T16Q2ForwardGroup<64, 51, 187, 34, 12, true>(
            tables, primary_base, recovery_base, recovery_count);
    }
}

template<unsigned Output, unsigned Log>
static LEO_FORCE_INLINE void T16Q4TailOutput(
    const unsigned char* tables,
    void* recovery_base,
    __m256i source0,
    __m256i source1)
{
    __m256i output0 = T16Q2Load(recovery_base, Output * 64U);
    __m256i output1 = T16Q2Load(recovery_base, Output * 64U + 32U);
    T16Q2MulAdd<Log>(tables, output0, source0);
    T16Q2MulAdd<Log>(tables, output1, source1);
    T16Q2Store(recovery_base, Output * 64U, output0);
    T16Q2Store(recovery_base, Output * 64U + 32U, output1);
}

static LEO_FORCE_INLINE void T16Q4AccumulateK65Tail(
    const unsigned char* tables,
    const void* data_base,
    void* recovery_base,
    unsigned recovery_count)
{
    LEO_DEBUG_ASSERT(recovery_count >= 9 && recovery_count <= 16);
    const __m256i source0 = T16Q2Load(data_base, 64U * 64U);
    const __m256i source1 = T16Q2Load(data_base, 64U * 64U + 32U);

    /* Exact [128,112] systematic Lagrange column for x=80 at y=0..15:

           L_x(y) = Z(y) / ((y + x) Z'(x)),
           Z(t) = product over s=16..127 of (t + s).

       Focused tests independently compare every output byte with the direct
       generator-matrix oracle rather than trusting these template logs. */
    T16Q4TailOutput<0, 124>(tables, recovery_base, source0, source1);
    T16Q4TailOutput<1, 248>(tables, recovery_base, source0, source1);
    T16Q4TailOutput<2, 72>(tables, recovery_base, source0, source1);
    T16Q4TailOutput<3, 227>(tables, recovery_base, source0, source1);
    T16Q4TailOutput<4, 143>(tables, recovery_base, source0, source1);
    T16Q4TailOutput<5, 199>(tables, recovery_base, source0, source1);
    T16Q4TailOutput<6, 62>(tables, recovery_base, source0, source1);
    T16Q4TailOutput<7, 132>(tables, recovery_base, source0, source1);
    T16Q4TailOutput<8, 33>(tables, recovery_base, source0, source1);
    if (recovery_count == 9)
        return;
    T16Q4TailOutput<9, 36>(tables, recovery_base, source0, source1);
    if (recovery_count == 10)
        return;
    T16Q4TailOutput<10, 31>(tables, recovery_base, source0, source1);
    if (recovery_count == 11)
        return;
    T16Q4TailOutput<11, 144>(tables, recovery_base, source0, source1);
    if (recovery_count == 12)
        return;
    T16Q4TailOutput<12, 66>(tables, recovery_base, source0, source1);
    if (recovery_count == 13)
        return;
    T16Q4TailOutput<13, 18>(tables, recovery_base, source0, source1);
    if (recovery_count == 14)
        return;
    T16Q4TailOutput<14, 9>(tables, recovery_base, source0, source1);
    if (recovery_count == 15)
        return;
    T16Q4TailOutput<15, 241>(tables, recovery_base, source0, source1);
}

#if defined(_MSC_VER)
#define LEO2_AVX2_T16_Q2_ENTRY __declspec(noinline)
#elif defined(__GNUC__) && !defined(__clang__) && defined(__ELF__)
#define LEO2_AVX2_T16_Q2_ENTRY \
    __attribute__((noinline, noipa, \
        section(".text.leo2_t16_q2"), aligned(64)))
#elif defined(__clang__) && defined(__ELF__)
#define LEO2_AVX2_T16_Q2_ENTRY \
    __attribute__((noinline, section(".text.leo2_t16_q2"), aligned(64)))
#elif defined(__GNUC__) || defined(__clang__)
#define LEO2_AVX2_T16_Q2_ENTRY __attribute__((noinline, aligned(64)))
#else
#define LEO2_AVX2_T16_Q2_ENTRY
#endif

/*
    Keep the larger Q3/Q4 circuit out of the ordinary .text wildcard used by
    the already-qualified Q2 kernels.  Our GNU-ld ELF release layout places
    this executable orphan after Leopard's earlier named hot sections, so the
    circuit does not displace mature encoder text merely by growing .text.
    Release-layout tests must verify the actual section order and addresses;
    other linkers receive only the noinline/alignment attributes below.
*/
#if defined(_MSC_VER)
#define LEO2_AVX2_T16_Q4_ENTRY __declspec(noinline)
#elif defined(__GNUC__) && !defined(__clang__) && defined(__ELF__)
#define LEO2_AVX2_T16_Q4_ENTRY \
    __attribute__((noinline, noipa, \
        section(".leo2_z_t16_q4_b64"), aligned(64)))
#elif defined(__clang__) && defined(__ELF__)
#define LEO2_AVX2_T16_Q4_ENTRY \
    __attribute__((noinline, \
        section(".leo2_z_t16_q4_b64"), aligned(64)))
#elif defined(__GNUC__) || defined(__clang__)
#define LEO2_AVX2_T16_Q4_ENTRY __attribute__((noinline, aligned(64)))
#else
#define LEO2_AVX2_T16_Q4_ENTRY
#endif

LEO2_AVX2_T16_Q2_ENTRY void AVX2FF8HighEncodeT16Q2B64Fused(
    const void* data_base,
    void* recovery_base,
    void* temporary_base,
    unsigned original_count,
    unsigned recovery_count)
{
    LEO_DEBUG_ASSERT(
        data_base != NULL && recovery_base != NULL && temporary_base != NULL);
    LEO_DEBUG_ASSERT(original_count >= 17 && original_count <= 32);
    LEO_DEBUG_ASSERT(recovery_count >= 9 && recovery_count <= 16);
    const unsigned char* tables =
        static_cast<const unsigned char*>(GetAVX2FF8Tables());
    LEO_DEBUG_ASSERT(tables != NULL);
    if (!tables)
        return;

    T16Q2EncodeSelected<64, 4>(tables, data_base, recovery_base,
        temporary_base, original_count, recovery_count);
}

LEO2_AVX2_T16_Q4_ENTRY void AVX2FF8HighEncodeT16Q4B64Fused(
    const void* data_base,
    void* recovery_base,
    void* temporary_base,
    unsigned original_count,
    unsigned recovery_count)
{
    LEO_DEBUG_ASSERT(
        data_base != NULL && recovery_base != NULL && temporary_base != NULL);
    LEO_DEBUG_ASSERT(original_count >= 33 && original_count <= 65);
    LEO_DEBUG_ASSERT(recovery_count >= 9 && recovery_count <= 16);
    const unsigned char* tables =
        static_cast<const unsigned char*>(GetAVX2FF8Tables());
    LEO_DEBUG_ASSERT(tables != NULL);
    if (!tables)
        return;

    T16Q4Encode(tables, data_base, recovery_base, temporary_base,
        original_count <= 64 ? original_count : 64U, recovery_count);
    if (original_count == 65)
    {
        T16Q4AccumulateK65Tail(
            tables, data_base, recovery_base, recovery_count);
    }
}

LEO2_AVX2_T16_Q2_ENTRY void AVX2FF8HighEncodeT16Q2B256Fused(
    const void* data_base,
    void* recovery_base,
    void* temporary_base,
    unsigned original_count,
    unsigned recovery_count)
{
    LEO_DEBUG_ASSERT(
        data_base != NULL && recovery_base != NULL && temporary_base != NULL);
    LEO_DEBUG_ASSERT(original_count >= 17 && original_count <= 32);
    LEO_DEBUG_ASSERT(recovery_count >= 9 && recovery_count <= 16);
    const unsigned char* tables =
        static_cast<const unsigned char*>(GetAVX2FF8Tables());
    LEO_DEBUG_ASSERT(tables != NULL);
    if (!tables)
        return;

    if (original_count <= 20)
        T16Q2EncodePartialSelected<256, 1>(
            tables, data_base, recovery_base,
            temporary_base, original_count, recovery_count);
    else if (original_count <= 24)
        T16Q2EncodePartialSelected<256, 2>(
            tables, data_base, recovery_base,
            temporary_base, original_count, recovery_count);
    else if (original_count <= 28)
        T16Q2EncodePartialSelected<256, 3>(
            tables, data_base, recovery_base,
            temporary_base, original_count, recovery_count);
    else
        T16Q2EncodeSelected<256, 4>(tables, data_base, recovery_base,
            temporary_base, original_count, recovery_count);
}

#undef LEO2_AVX2_T16_Q2_ENTRY
#undef LEO2_AVX2_T16_Q4_ENTRY

#endif

}} // namespace leopard::backend
