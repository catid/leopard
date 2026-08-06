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
static const uint64_t kT16Q2ShardBytes = 64;

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
        for (uint64_t offset = 0; offset < kT16Q2ShardBytes; offset += 32)
        {
            T16Q2Store(work_base,
                (OutputRow + 0U) * kT16Q2ShardBytes + offset, zero);
            T16Q2Store(work_base,
                (OutputRow + 1U) * kT16Q2ShardBytes + offset, zero);
            T16Q2Store(work_base,
                (OutputRow + 2U) * kT16Q2ShardBytes + offset, zero);
            T16Q2Store(work_base,
                (OutputRow + 3U) * kT16Q2ShardBytes + offset, zero);
        }
        return;
    }

    for (uint64_t offset = 0; offset < kT16Q2ShardBytes; offset += 32)
    {
        const __m256i zero = _mm256_setzero_si256();
        __m256i value0 = !PartialInput || InputRow + 0U < original_count
            ? T16Q2Load(data_base,
                (InputRow + 0U) * kT16Q2ShardBytes + offset)
            : zero;
        __m256i value1 = !PartialInput || InputRow + 1U < original_count
            ? T16Q2Load(data_base,
                (InputRow + 1U) * kT16Q2ShardBytes + offset)
            : zero;
        __m256i value2 = !PartialInput || InputRow + 2U < original_count
            ? T16Q2Load(data_base,
                (InputRow + 2U) * kT16Q2ShardBytes + offset)
            : zero;
        __m256i value3 = !PartialInput || InputRow + 3U < original_count
            ? T16Q2Load(data_base,
                (InputRow + 3U) * kT16Q2ShardBytes + offset)
            : zero;
        T16Q2IFFT4<Log01, Log23, Log02>(
            tables, value0, value1, value2, value3);
        T16Q2Store(work_base,
            (OutputRow + 0U) * kT16Q2ShardBytes + offset, value0);
        T16Q2Store(work_base,
            (OutputRow + 1U) * kT16Q2ShardBytes + offset, value1);
        T16Q2Store(work_base,
            (OutputRow + 2U) * kT16Q2ShardBytes + offset, value2);
        T16Q2Store(work_base,
            (OutputRow + 3U) * kT16Q2ShardBytes + offset, value3);
    }
}

static LEO_FORCE_INLINE void T16Q2FusedOuter(
    const unsigned char* tables,
    void* recovery_base,
    void* temporary_base)
{
    for (unsigned column = 0; column < 4; ++column)
    {
        for (uint64_t offset = 0; offset < kT16Q2ShardBytes; offset += 32)
        {
            __m256i first0 = T16Q2Load(
                recovery_base,
                (column + 0U) * kT16Q2ShardBytes + offset);
            __m256i first1 = T16Q2Load(
                recovery_base,
                (column + 4U) * kT16Q2ShardBytes + offset);
            __m256i first2 = T16Q2Load(
                recovery_base,
                (column + 8U) * kT16Q2ShardBytes + offset);
            __m256i first3 = T16Q2Load(
                recovery_base,
                (column + 12U) * kT16Q2ShardBytes + offset);
            __m256i second0 = T16Q2Load(
                temporary_base,
                (column + 0U) * kT16Q2ShardBytes + offset);
            __m256i second1 = T16Q2Load(
                temporary_base,
                (column + 4U) * kT16Q2ShardBytes + offset);
            __m256i second2 = T16Q2Load(
                temporary_base,
                (column + 8U) * kT16Q2ShardBytes + offset);
            __m256i second3 = T16Q2Load(
                temporary_base,
                (column + 12U) * kT16Q2ShardBytes + offset);

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
                (column + 0U) * kT16Q2ShardBytes + offset, first0);
            T16Q2Store(
                recovery_base,
                (column + 4U) * kT16Q2ShardBytes + offset, first1);
            T16Q2Store(
                recovery_base,
                (column + 8U) * kT16Q2ShardBytes + offset, first2);
            T16Q2Store(
                recovery_base,
                (column + 12U) * kT16Q2ShardBytes + offset, first3);
        }
    }
}

template<
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
    for (uint64_t offset = 0; offset < kT16Q2ShardBytes; offset += 32)
    {
        __m256i value0 = T16Q2Load(
            work_base,
            (BaseRow + 0U) * kT16Q2ShardBytes + offset);
        __m256i value1 = T16Q2Load(
            work_base,
            (BaseRow + 1U) * kT16Q2ShardBytes + offset);
        __m256i value2 = T16Q2Load(
            work_base,
            (BaseRow + 2U) * kT16Q2ShardBytes + offset);
        __m256i value3 = T16Q2Load(
            work_base,
            (BaseRow + 3U) * kT16Q2ShardBytes + offset);
        T16Q2FFT4<Log01, Log23, Log02>(
            tables, value0, value1, value2, value3);
        if (!PartialOutput || BaseRow + 0U < recovery_count)
            T16Q2Store(recovery_base,
                (BaseRow + 0U) * kT16Q2ShardBytes + offset, value0);
        if (!PartialOutput || BaseRow + 1U < recovery_count)
            T16Q2Store(recovery_base,
                (BaseRow + 1U) * kT16Q2ShardBytes + offset, value1);
        if (!PartialOutput || BaseRow + 2U < recovery_count)
            T16Q2Store(recovery_base,
                (BaseRow + 2U) * kT16Q2ShardBytes + offset, value2);
        if (!PartialOutput || BaseRow + 3U < recovery_count)
            T16Q2Store(recovery_base,
                (BaseRow + 3U) * kT16Q2ShardBytes + offset, value3);
    }
}

template<bool PartialInput, bool PartialOutput>
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
        ? temporary_bytes + 16U * kT16Q2ShardBytes : temporary_base;

    /* First two inverse layers, FFTSkewStorage + 16. */
    T16Q2InverseGroup<219, 7, 153, 0, 0, false>(
        tables, data_base, primary_base, original_count);
    T16Q2InverseGroup<111, 28, 102, 4, 4, false>(
        tables, data_base, primary_base, original_count);
    T16Q2InverseGroup<183, 224, 51, 8, 8, false>(
        tables, data_base, primary_base, original_count);
    T16Q2InverseGroup<131, 222, 187, 12, 12, false>(
        tables, data_base, primary_base, original_count);

    /* First two inverse layers, FFTSkewStorage + 32. */
    T16Q2InverseGroup<196, 76, 219, 16, 0, PartialInput>(
        tables, data_base, secondary_base, original_count);
    T16Q2InverseGroup<54, 99, 7, 20, 4, PartialInput>(
        tables, data_base, secondary_base, original_count);
    T16Q2InverseGroup<19, 49, 111, 24, 8, PartialInput>(
        tables, data_base, secondary_base, original_count);
    T16Q2InverseGroup<67, 52, 28, 28, 12, PartialInput>(
        tables, data_base, secondary_base, original_count);

    T16Q2FusedOuter(tables, primary_base, secondary_base);

    /* Final forward layer, FFTSkewStorage + 0. */
    T16Q2ForwardGroup<255, 85, 255, 0, PartialOutput>(
        tables, primary_base, recovery_base, recovery_count);
    T16Q2ForwardGroup<17, 34, 85, 4, PartialOutput>(
        tables, primary_base, recovery_base, recovery_count);
    T16Q2ForwardGroup<153, 102, 17, 8, PartialOutput>(
        tables, primary_base, recovery_base, recovery_count);
    T16Q2ForwardGroup<51, 187, 34, 12, PartialOutput>(
        tables, primary_base, recovery_base, recovery_count);
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

    if (original_count == 32 && recovery_count == 16)
        T16Q2Encode<false, false>(tables, data_base, recovery_base,
            temporary_base, original_count, recovery_count);
    else if (recovery_count == 16)
        T16Q2Encode<true, false>(tables, data_base, recovery_base,
            temporary_base, original_count, recovery_count);
    else
        T16Q2Encode<true, true>(tables, data_base, recovery_base,
            temporary_base, original_count, recovery_count);
}

#undef LEO2_AVX2_T16_Q2_ENTRY

#endif

}} // namespace leopard::backend
