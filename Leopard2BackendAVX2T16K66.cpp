/*
    Copyright (c) 2017 Christopher A. Taylor.  All rights reserved.
    See Leopard2Backend.h for the full BSD license notice.
*/

#include "Leopard2Backend.h"
#include "LeopardCommon.h"

#include <immintrin.h>

namespace leopard { namespace backend {

#if LEO2_EXPERIMENT_HIGH_T16_Q2_B64_FUSED && defined(LEO_HAS_FF8) && \
    !defined(LEO2_AVX512_VARIANT) && !defined(LEO2_GFNI_VARIANT)

namespace {

static LEO_FORCE_INLINE __m256i Broadcast(const uint8_t* table)
{
    return _mm256_broadcastsi128_si256(_mm_loadu_si128(
        reinterpret_cast<const __m128i*>(table)));
}

static LEO_FORCE_INLINE __m256i Product(
    __m256i data,
    __m256i low_table,
    __m256i high_table)
{
    const __m256i nibble_mask = _mm256_set1_epi8(15);
    const __m256i low = _mm256_shuffle_epi8(
        low_table, _mm256_and_si256(data, nibble_mask));
    const __m256i high = _mm256_shuffle_epi8(
        high_table,
        _mm256_and_si256(_mm256_srli_epi64(data, 4), nibble_mask));
    return _mm256_xor_si256(low, high);
}

template<unsigned Output, unsigned Log>
static LEO_FORCE_INLINE void AccumulateTailOutput(
    const uint8_t* tables,
    const void* data_base,
    void* recovery_base)
{
    const uint8_t* const table = tables + Log * 32U;
    const __m256i low_table = Broadcast(table);
    const __m256i high_table = Broadcast(table + 16U);
    for (unsigned offset = 0; offset < 64; offset += 32)
    {
        const __m256i source = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(
                static_cast<const uint8_t*>(data_base) + 65U * 64U +
                offset));
        __m256i output = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(
                static_cast<uint8_t*>(recovery_base) + Output * 64U +
                offset));
        output = _mm256_xor_si256(
            output, Product(source, low_table, high_table));
        _mm256_storeu_si256(
            reinterpret_cast<__m256i*>(
                static_cast<uint8_t*>(recovery_base) + Output * 64U +
                offset), output);
    }
}

} // namespace

#if defined(_MSC_VER)
#define LEO2_AVX2_T16_K66_ENTRY __declspec(noinline)
#elif defined(__GNUC__) && !defined(__clang__) && defined(__ELF__)
#define LEO2_AVX2_T16_K66_ENTRY \
    __attribute__((noinline, noipa, \
        section(".leo2_z_t16_k66_b64"), aligned(64)))
#elif defined(__clang__) && defined(__ELF__)
#define LEO2_AVX2_T16_K66_ENTRY \
    __attribute__((noinline, \
        section(".leo2_z_t16_k66_b64"), aligned(64)))
#elif defined(__GNUC__) || defined(__clang__)
#define LEO2_AVX2_T16_K66_ENTRY __attribute__((noinline, aligned(64)))
#else
#define LEO2_AVX2_T16_K66_ENTRY
#endif

LEO2_AVX2_T16_K66_ENTRY void AVX2FF8HighEncodeT16K66R16B64Fused(
    const void* data_base,
    void* recovery_base,
    void* temporary_base)
{
    LEO_DEBUG_ASSERT(
        data_base != NULL && recovery_base != NULL && temporary_base != NULL);
    const uint8_t* const tables =
        static_cast<const uint8_t*>(GetAVX2FF8Tables());
    LEO_DEBUG_ASSERT(tables != NULL);
    if (!tables)
        return;

    /* Preserve the already-qualified K=33..65 Q4 entry byte-for-byte.  K66
       is exactly its K65 result plus the systematic Lagrange column at
       coordinate 81 in the same [128,112] parent. */
    AVX2FF8HighEncodeT16Q4B64Fused(
        data_base, recovery_base, temporary_base, 65, 16);

    /* Two independent derivations reproduce the existing coordinate-80
       column and agree on this coordinate-81 column.  Its logs are the
       expected adjacent-pair permutation of the K65 logs. */
    AccumulateTailOutput<0, 248>(tables, data_base, recovery_base);
    AccumulateTailOutput<1, 124>(tables, data_base, recovery_base);
    AccumulateTailOutput<2, 227>(tables, data_base, recovery_base);
    AccumulateTailOutput<3, 72>(tables, data_base, recovery_base);
    AccumulateTailOutput<4, 199>(tables, data_base, recovery_base);
    AccumulateTailOutput<5, 143>(tables, data_base, recovery_base);
    AccumulateTailOutput<6, 132>(tables, data_base, recovery_base);
    AccumulateTailOutput<7, 62>(tables, data_base, recovery_base);
    AccumulateTailOutput<8, 36>(tables, data_base, recovery_base);
    AccumulateTailOutput<9, 33>(tables, data_base, recovery_base);
    AccumulateTailOutput<10, 144>(tables, data_base, recovery_base);
    AccumulateTailOutput<11, 31>(tables, data_base, recovery_base);
    AccumulateTailOutput<12, 18>(tables, data_base, recovery_base);
    AccumulateTailOutput<13, 66>(tables, data_base, recovery_base);
    AccumulateTailOutput<14, 241>(tables, data_base, recovery_base);
    AccumulateTailOutput<15, 9>(tables, data_base, recovery_base);
}

#undef LEO2_AVX2_T16_K66_ENTRY

#endif

}} // namespace leopard::backend
