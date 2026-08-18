/*
    Copyright (c) 2017 Christopher A. Taylor.  All rights reserved.
    See Leopard2Backend.h for the full BSD license notice.
*/

#include "Leopard2Backend.h"
#include "LeopardCommon.h"

#include <immintrin.h>

namespace leopard { namespace backend {

#if LEO2_EXPERIMENT_HIGH_T8_TWO_BLOCK_BINDING && \
    defined(LEO_HAS_FF8) && \
    !defined(LEO2_DIAGNOSTIC_DISABLE_HIGH_T8_VECTOR) && \
    !defined(LEO2_AVX512_VARIANT) && !defined(LEO2_GFNI_VARIANT)

/* Published by the ordinary AVX2 translation unit after table setup. */
const void* GetAVX2FF8Tables();

namespace {

static const unsigned kNibbleTableBytes = 32U;

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

static LEO_FORCE_INLINE const uint8_t* Table(
    const uint8_t* tables,
    uint16_t log)
{
    return tables + static_cast<unsigned>(log) * kNibbleTableBytes;
}

static LEO_FORCE_INLINE void MulAddPrepared(
    __m256i& destination,
    __m256i source,
    __m256i low_table,
    __m256i high_table)
{
    destination = _mm256_xor_si256(
        destination, Product(source, low_table, high_table));
}

static LEO_FORCE_INLINE void MulAdd(
    const uint8_t* tables,
    __m256i& destination,
    __m256i source,
    uint16_t log)
{
    const uint8_t* const table = Table(tables, log);
    MulAddPrepared(
        destination, source, Broadcast(table), Broadcast(table + 16));
}

static LEO_FORCE_INLINE void IFFTRadix4(
    const uint8_t* tables,
    __m256i& value0,
    __m256i& value1,
    __m256i& value2,
    __m256i& value3,
    uint16_t log01,
    uint16_t log23,
    uint16_t log02)
{
    static const uint16_t kZeroSkew = 255;
    value1 = _mm256_xor_si256(value1, value0);
    if (log01 != kZeroSkew)
        MulAdd(tables, value0, value1, log01);
    value3 = _mm256_xor_si256(value3, value2);
    if (log23 != kZeroSkew)
        MulAdd(tables, value2, value3, log23);
    value2 = _mm256_xor_si256(value2, value0);
    value3 = _mm256_xor_si256(value3, value1);
    if (log02 != kZeroSkew)
    {
        const uint8_t* const table = Table(tables, log02);
        const __m256i low_table = Broadcast(table);
        const __m256i high_table = Broadcast(table + 16);
        MulAddPrepared(value0, value2, low_table, high_table);
        MulAddPrepared(value1, value3, low_table, high_table);
    }
}

static LEO_FORCE_INLINE void IFFTDistance4(
    const uint8_t* tables,
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
    value4 = _mm256_xor_si256(value4, value0);
    value5 = _mm256_xor_si256(value5, value1);
    value6 = _mm256_xor_si256(value6, value2);
    value7 = _mm256_xor_si256(value7, value3);
    if (log == kZeroSkew)
        return;
    const uint8_t* const table = Table(tables, log);
    const __m256i low_table = Broadcast(table);
    const __m256i high_table = Broadcast(table + 16);
    MulAddPrepared(value0, value4, low_table, high_table);
    MulAddPrepared(value1, value5, low_table, high_table);
    MulAddPrepared(value2, value6, low_table, high_table);
    MulAddPrepared(value3, value7, low_table, high_table);
}

static LEO_FORCE_INLINE void InverseValues(
    const uint8_t* tables,
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
    IFFTRadix4(tables, value0, value1, value2, value3,
        inverse_skew[1], inverse_skew[3], inverse_skew[2]);
    IFFTRadix4(tables, value4, value5, value6, value7,
        inverse_skew[5], inverse_skew[7], inverse_skew[6]);
    IFFTDistance4(tables,
        value0, value1, value2, value3,
        value4, value5, value6, value7, inverse_skew[4]);
}

static LEO_FORCE_INLINE void FFTRadix4Distance2(
    const uint8_t* tables,
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
        const uint8_t* const table = Table(tables, log02);
        const __m256i low_table = Broadcast(table);
        const __m256i high_table = Broadcast(table + 16);
        MulAddPrepared(value0, value4, low_table, high_table);
        MulAddPrepared(value1, value5, low_table, high_table);
        MulAddPrepared(value2, value6, low_table, high_table);
        MulAddPrepared(value3, value7, low_table, high_table);
    }
    value4 = _mm256_xor_si256(value4, value0);
    value5 = _mm256_xor_si256(value5, value1);
    value6 = _mm256_xor_si256(value6, value2);
    value7 = _mm256_xor_si256(value7, value3);
    if (log01 != kZeroSkew)
    {
        const uint8_t* const table = Table(tables, log01);
        const __m256i low_table = Broadcast(table);
        const __m256i high_table = Broadcast(table + 16);
        MulAddPrepared(value0, value2, low_table, high_table);
        MulAddPrepared(value1, value3, low_table, high_table);
    }
    value2 = _mm256_xor_si256(value2, value0);
    value3 = _mm256_xor_si256(value3, value1);
    if (log23 != kZeroSkew)
    {
        const uint8_t* const table = Table(tables, log23);
        const __m256i low_table = Broadcast(table);
        const __m256i high_table = Broadcast(table + 16);
        MulAddPrepared(value4, value6, low_table, high_table);
        MulAddPrepared(value5, value7, low_table, high_table);
    }
    value6 = _mm256_xor_si256(value6, value4);
    value7 = _mm256_xor_si256(value7, value5);
}

static LEO_FORCE_INLINE void FFT2(
    const uint8_t* tables,
    __m256i& value0,
    __m256i& value1,
    uint16_t log)
{
    if (log != 255)
        MulAdd(tables, value0, value1, log);
    value1 = _mm256_xor_si256(value1, value0);
}

static LEO_FORCE_INLINE void ForwardValues(
    const uint8_t* tables,
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
    FFTRadix4Distance2(tables,
        value0, value1, value2, value3,
        value4, value5, value6, value7,
        forward_skew[2], forward_skew[6], forward_skew[4]);
    FFT2(tables, value0, value1, forward_skew[1]);
    FFT2(tables, value2, value3, forward_skew[3]);
    FFT2(tables, value4, value5, forward_skew[5]);
    FFT2(tables, value6, value7, forward_skew[7]);
}

template<unsigned BlockBase, unsigned LiveRows>
static LEO_FORCE_INLINE void LoadInversePacked(
    const uint8_t* tables,
    const void* const* data,
    uint64_t offset,
    const uint8_t* inverse_skew,
    __m256i& value0,
    __m256i& value1,
    __m256i& value2,
    __m256i& value3,
    __m256i& value4,
    __m256i& value5,
    __m256i& value6,
    __m256i& value7)
{
    static_assert(LiveRows == 6 || LiveRows == 8,
        "K62 T8 blocks are complete or have six live rows");
#define LEO2_T8_K62_LOAD(Lane) \
    _mm256_loadu_si256(reinterpret_cast<const __m256i*>( \
        static_cast<const uint8_t*>(data[BlockBase + Lane]) + offset))
    value0 = LEO2_T8_K62_LOAD(0);
    value1 = LEO2_T8_K62_LOAD(1);
    value2 = LEO2_T8_K62_LOAD(2);
    value3 = LEO2_T8_K62_LOAD(3);
    value4 = LEO2_T8_K62_LOAD(4);
    value5 = LEO2_T8_K62_LOAD(5);
    if (LiveRows == 8)
    {
        value6 = LEO2_T8_K62_LOAD(6);
        value7 = LEO2_T8_K62_LOAD(7);
    }
    else
    {
        value6 = _mm256_setzero_si256();
        value7 = _mm256_setzero_si256();
    }
#undef LEO2_T8_K62_LOAD
    InverseValues(tables,
        value0, value1, value2, value3,
        value4, value5, value6, value7, inverse_skew);
}

template<bool Accumulate, bool FinishForward>
static LEO_FORCE_INLINE void StorePacked(
    const uint8_t* tables,
    void* const* recovery,
    uint64_t offset,
    const uint8_t* forward_skew,
    __m256i value0,
    __m256i value1,
    __m256i value2,
    __m256i value3,
    __m256i value4,
    __m256i value5,
    __m256i value6,
    __m256i value7)
{
#define LEO2_T8_K62_OUTPUT(Lane) \
    reinterpret_cast<__m256i*>( \
        static_cast<uint8_t*>(recovery[Lane]) + offset)
    if (Accumulate)
    {
        value0 = _mm256_xor_si256(
            value0, _mm256_loadu_si256(LEO2_T8_K62_OUTPUT(0)));
        value1 = _mm256_xor_si256(
            value1, _mm256_loadu_si256(LEO2_T8_K62_OUTPUT(1)));
        value2 = _mm256_xor_si256(
            value2, _mm256_loadu_si256(LEO2_T8_K62_OUTPUT(2)));
        value3 = _mm256_xor_si256(
            value3, _mm256_loadu_si256(LEO2_T8_K62_OUTPUT(3)));
        value4 = _mm256_xor_si256(
            value4, _mm256_loadu_si256(LEO2_T8_K62_OUTPUT(4)));
        value5 = _mm256_xor_si256(
            value5, _mm256_loadu_si256(LEO2_T8_K62_OUTPUT(5)));
        value6 = _mm256_xor_si256(
            value6, _mm256_loadu_si256(LEO2_T8_K62_OUTPUT(6)));
        value7 = _mm256_xor_si256(
            value7, _mm256_loadu_si256(LEO2_T8_K62_OUTPUT(7)));
    }
    if (FinishForward)
    {
        ForwardValues(tables,
            value0, value1, value2, value3,
            value4, value5, value6, value7, forward_skew);
    }
    _mm256_storeu_si256(LEO2_T8_K62_OUTPUT(0), value0);
    _mm256_storeu_si256(LEO2_T8_K62_OUTPUT(1), value1);
    _mm256_storeu_si256(LEO2_T8_K62_OUTPUT(2), value2);
    _mm256_storeu_si256(LEO2_T8_K62_OUTPUT(3), value3);
    _mm256_storeu_si256(LEO2_T8_K62_OUTPUT(4), value4);
    _mm256_storeu_si256(LEO2_T8_K62_OUTPUT(5), value5);
    _mm256_storeu_si256(LEO2_T8_K62_OUTPUT(6), value6);
    _mm256_storeu_si256(LEO2_T8_K62_OUTPUT(7), value7);
#undef LEO2_T8_K62_OUTPUT
}

} // namespace

#if defined(_MSC_VER)
#define LEO2_AVX2_T8_K62_B64_ENTRY __declspec(noinline)
#elif defined(__GNUC__) && !defined(__clang__) && defined(__ELF__)
#define LEO2_AVX2_T8_K62_B64_ENTRY \
    __attribute__((noinline, noipa, \
        section(".leo2_z_t8_k62_b64_fused"), aligned(64)))
#elif defined(__clang__) && defined(__ELF__)
#define LEO2_AVX2_T8_K62_B64_ENTRY \
    __attribute__((noinline, \
        section(".leo2_z_t8_k62_b64_fused"), aligned(64)))
#elif defined(__GNUC__) || defined(__clang__)
#define LEO2_AVX2_T8_K62_B64_ENTRY __attribute__((noinline, aligned(64)))
#else
#define LEO2_AVX2_T8_K62_B64_ENTRY
#endif

LEO2_AVX2_T8_K62_B64_ENTRY bool AVX2FF8HighEncodeK62R8T8B64Fused(
    const void* const* data,
    void* const* recovery,
    const uint8_t* first_inverse_skew,
    const uint8_t* forward_skew)
{
    LEO_DEBUG_ASSERT(
        data != NULL && recovery != NULL &&
        first_inverse_skew != NULL && forward_skew != NULL);
    const uint8_t* const tables =
        static_cast<const uint8_t*>(GetAVX2FF8Tables());
    LEO_DEBUG_ASSERT(tables != NULL);
    if (!tables)
        return false;
    for (uint64_t offset = 0; offset < 64; offset += 32)
    {
        __m256i value0;
        __m256i value1;
        __m256i value2;
        __m256i value3;
        __m256i value4;
        __m256i value5;
        __m256i value6;
        __m256i value7;
#define LEO2_T8_K62_BLOCK(Base, Live, Accumulate, Finish) \
        LoadInversePacked<Base, Live>( \
            tables, data, offset, first_inverse_skew + Base, \
            value0, value1, value2, value3, \
            value4, value5, value6, value7); \
        StorePacked<Accumulate, Finish>( \
            tables, recovery, offset, forward_skew, \
            value0, value1, value2, value3, \
            value4, value5, value6, value7)
        LEO2_T8_K62_BLOCK(0, 8, false, false);
        LEO2_T8_K62_BLOCK(8, 8, true, false);
        LEO2_T8_K62_BLOCK(16, 8, true, false);
        LEO2_T8_K62_BLOCK(24, 8, true, false);
        LEO2_T8_K62_BLOCK(32, 8, true, false);
        LEO2_T8_K62_BLOCK(40, 8, true, false);
        LEO2_T8_K62_BLOCK(48, 8, true, false);
        LEO2_T8_K62_BLOCK(56, 6, true, true);
#undef LEO2_T8_K62_BLOCK
    }
    return true;
}

#undef LEO2_AVX2_T8_K62_B64_ENTRY

#endif

}} // namespace leopard::backend
