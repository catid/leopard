/*
    Copyright (c) 2017 Christopher A. Taylor.  All rights reserved.
    See Leopard2Backend.h for the full BSD license notice.
*/

#include "Leopard2Backend.h"
#include "LeopardCommon.h"

#include <immintrin.h>

namespace leopard { namespace backend {

#if defined(LEO2_EXPERIMENT_HIGH_T32_B256_GENERATED) && \
    defined(LEO2_HAVE_AVX2_BACKEND) && !defined(NO_LEO_HAS_FF8)

namespace {

#ifndef LEO2_DIAGNOSTIC_DISABLE_HIGH_T32_B256_GENERATED
#define LEO2_DIAGNOSTIC_DISABLE_HIGH_T32_B256_GENERATED 0
#endif

/*
    A volatile data selector keeps the control and candidate instruction
    streams identical.  Diagnostic builds change this initializer only; when
    disabled, the public routing TU falls through to the mature transform.
*/
static volatile uint32_t g_t32_b256_generated_mode =
    1U + LEO2_DIAGNOSTIC_DISABLE_HIGH_T32_B256_GENERATED;

/*
    The ordinary TU stores two consecutive 16-byte rows per logarithm.  Read
    its object representation through bytes instead of aliasing it as a second
    private struct type across translation units.
*/
static LEO_FORCE_INLINE const uint8_t* Tables()
{
    return static_cast<const uint8_t*>(GetAVX2FF8Tables());
}

static LEO_FORCE_INLINE __m256i Broadcast(const uint8_t table[16])
{
    return _mm256_broadcastsi128_si256(_mm_loadu_si128(
        reinterpret_cast<const __m128i*>(table)));
}

static LEO_FORCE_INLINE __m256i Product(
    __m256i data,
    __m256i low_table,
    __m256i high_table)
{
    const __m256i mask = _mm256_set1_epi8(15);
    const __m256i low = _mm256_shuffle_epi8(
        low_table, _mm256_and_si256(data, mask));
    const __m256i high = _mm256_shuffle_epi8(high_table,
        _mm256_and_si256(_mm256_srli_epi64(data, 4), mask));
    return _mm256_xor_si256(low, high);
}

static LEO_FORCE_INLINE __m256i Load(const void* pointer, unsigned offset)
{
    const uint8_t* const bytes = static_cast<const uint8_t*>(pointer) + offset;
    return _mm256_loadu_si256(reinterpret_cast<const __m256i*>(bytes));
}

static LEO_FORCE_INLINE void Store(
    void* pointer,
    unsigned offset,
    __m256i value)
{
    uint8_t* const bytes = static_cast<uint8_t*>(pointer) + offset;
    _mm256_storeu_si256(reinterpret_cast<__m256i*>(bytes), value);
}

template<unsigned Log>
static LEO_FORCE_INLINE void MulAdd(
    const uint8_t* tables,
    __m256i& destination,
    __m256i source)
{
    const uint8_t* const table = tables + Log * 32U;
    destination = _mm256_xor_si256(destination,
        Product(source, Broadcast(table), Broadcast(table + 16)));
}

template<unsigned Log>
static LEO_FORCE_INLINE void IFFT2(
    const uint8_t* tables,
    __m256i& x,
    __m256i& y)
{
    y = _mm256_xor_si256(y, x);
    if (Log != 255)
        MulAdd<Log>(tables, x, y);
}

template<unsigned Log>
static LEO_FORCE_INLINE void FFT2(
    const uint8_t* tables,
    __m256i& x,
    __m256i& y)
{
    if (Log != 255)
        MulAdd<Log>(tables, x, y);
    y = _mm256_xor_si256(y, x);
}

template<unsigned Log01, unsigned Log23, unsigned Log02>
static LEO_FORCE_INLINE void IFFT4(
    const uint8_t* tables,
    __m256i& x0, __m256i& x1, __m256i& x2, __m256i& x3)
{
    IFFT2<Log01>(tables, x0, x1);
    IFFT2<Log23>(tables, x2, x3);
    IFFT2<Log02>(tables, x0, x2);
    IFFT2<Log02>(tables, x1, x3);
}

template<unsigned Log01, unsigned Log23, unsigned Log02>
static LEO_FORCE_INLINE void FFT4(
    const uint8_t* tables,
    __m256i& x0, __m256i& x1, __m256i& x2, __m256i& x3)
{
    FFT2<Log02>(tables, x0, x2);
    FFT2<Log02>(tables, x1, x3);
    FFT2<Log01>(tables, x0, x1);
    FFT2<Log23>(tables, x2, x3);
}

template<unsigned Log01, unsigned Log23, unsigned Log45,
    unsigned Log67, unsigned Log02, unsigned Log46, unsigned Log04>
static LEO_FORCE_INLINE void InverseGroup(
    const uint8_t* tables,
    const void* const* data,
    void* const* work,
    unsigned base,
    unsigned offset)
{
    __m256i x0 = Load(data[base], offset);
    __m256i x1 = Load(data[base + 1], offset);
    __m256i x2 = Load(data[base + 2], offset);
    __m256i x3 = Load(data[base + 3], offset);
    __m256i x4 = Load(data[base + 4], offset);
    __m256i x5 = Load(data[base + 5], offset);
    __m256i x6 = Load(data[base + 6], offset);
    __m256i x7 = Load(data[base + 7], offset);
    IFFT2<Log01>(tables, x0, x1);
    IFFT2<Log23>(tables, x2, x3);
    IFFT2<Log45>(tables, x4, x5);
    IFFT2<Log67>(tables, x6, x7);
    IFFT2<Log02>(tables, x0, x2);
    IFFT2<Log02>(tables, x1, x3);
    IFFT2<Log46>(tables, x4, x6);
    IFFT2<Log46>(tables, x5, x7);
    IFFT2<Log04>(tables, x0, x4);
    IFFT2<Log04>(tables, x1, x5);
    IFFT2<Log04>(tables, x2, x6);
    IFFT2<Log04>(tables, x3, x7);
    Store(work[base], offset, x0);
    Store(work[base + 1], offset, x1);
    Store(work[base + 2], offset, x2);
    Store(work[base + 3], offset, x3);
    Store(work[base + 4], offset, x4);
    Store(work[base + 5], offset, x5);
    Store(work[base + 6], offset, x6);
    Store(work[base + 7], offset, x7);
}

static LEO_FORCE_INLINE void OuterGroups(
    const uint8_t* tables,
    void* const* work,
    unsigned column,
    unsigned offset)
{
    __m256i x0 = Load(work[column], offset);
    __m256i x1 = Load(work[column + 8], offset);
    __m256i x2 = Load(work[column + 16], offset);
    __m256i x3 = Load(work[column + 24], offset);
    IFFT4<17, 34, 85>(tables, x0, x1, x2, x3);
    FFT4<255, 85, 255>(tables, x0, x1, x2, x3);
    Store(work[column], offset, x0);
    Store(work[column + 8], offset, x1);
    Store(work[column + 16], offset, x2);
    Store(work[column + 24], offset, x3);
}

template<unsigned Log01, unsigned Log23, unsigned Log45,
    unsigned Log67, unsigned Log02, unsigned Log46, unsigned Log04>
static LEO_FORCE_INLINE void ForwardGroup(
    const uint8_t* tables,
    void* const* work,
    unsigned base,
    unsigned offset)
{
    __m256i x0 = Load(work[base], offset);
    __m256i x1 = Load(work[base + 1], offset);
    __m256i x2 = Load(work[base + 2], offset);
    __m256i x3 = Load(work[base + 3], offset);
    __m256i x4 = Load(work[base + 4], offset);
    __m256i x5 = Load(work[base + 5], offset);
    __m256i x6 = Load(work[base + 6], offset);
    __m256i x7 = Load(work[base + 7], offset);
    FFT2<Log04>(tables, x0, x4);
    FFT2<Log04>(tables, x1, x5);
    FFT2<Log04>(tables, x2, x6);
    FFT2<Log04>(tables, x3, x7);
    FFT2<Log02>(tables, x0, x2);
    FFT2<Log02>(tables, x1, x3);
    FFT2<Log46>(tables, x4, x6);
    FFT2<Log46>(tables, x5, x7);
    FFT2<Log01>(tables, x0, x1);
    FFT2<Log23>(tables, x2, x3);
    FFT2<Log45>(tables, x4, x5);
    FFT2<Log67>(tables, x6, x7);
    Store(work[base], offset, x0);
    Store(work[base + 1], offset, x1);
    Store(work[base + 2], offset, x2);
    Store(work[base + 3], offset, x3);
    Store(work[base + 4], offset, x4);
    Store(work[base + 5], offset, x5);
    Store(work[base + 6], offset, x6);
    Store(work[base + 7], offset, x7);
}

} // namespace

#if defined(_MSC_VER)
#define LEO2_AVX2_T32_B256_ENTRY __declspec(noinline)
#elif defined(__GNUC__) && !defined(__clang__) && defined(__ELF__)
#define LEO2_AVX2_T32_B256_ENTRY \
    __attribute__((noinline, noipa, section(".text.leo2_t32_b256"), \
        aligned(64)))
#elif defined(__clang__) && defined(__ELF__)
#define LEO2_AVX2_T32_B256_ENTRY \
    __attribute__((noinline, section(".text.leo2_t32_b256"), aligned(64)))
#elif defined(__GNUC__) || defined(__clang__)
#define LEO2_AVX2_T32_B256_ENTRY __attribute__((noinline, aligned(64)))
#else
#define LEO2_AVX2_T32_B256_ENTRY
#endif

LEO2_AVX2_T32_B256_ENTRY bool TryAVX2FF8HighEncodeT32B256(
    const void* const* data,
    void* const* recovery)
{
    if (g_t32_b256_generated_mode != 1U)
        return false;

    LEO_DEBUG_ASSERT(data != NULL && recovery != NULL);
    const uint8_t* const tables = Tables();
    LEO_DEBUG_ASSERT(tables != NULL);
    if (!tables)
        return false;

    /* Pass 1: inverse distances 1, 2, and 4. */
    for (unsigned offset = 0; offset < 256; offset += 32)
    {
        InverseGroup<196,76,54,99,219,7,153>(
            tables, data, recovery, 0, offset);
        InverseGroup<19,49,67,52,111,28,102>(
            tables, data, recovery, 8, offset);
        InverseGroup<137,226,26,108,183,224,51>(
            tables, data, recovery, 16, offset);
        InverseGroup<27,134,38,139,131,222,187>(
            tables, data, recovery, 24, offset);
    }

    /* Pass 2: inverse distances 8/16, then forward distances 16/8. */
    for (unsigned offset = 0; offset < 256; offset += 32)
        for (unsigned column = 0; column < 8; ++column)
            OuterGroups(tables, recovery, column, offset);

    /* Pass 3: forward distances 4, 2, and 1. */
    for (unsigned offset = 0; offset < 256; offset += 32)
    {
        ForwardGroup<255,85,17,34,255,85,255>(
            tables, recovery, 0, offset);
        ForwardGroup<153,102,51,187,17,34,85>(
            tables, recovery, 8, offset);
        ForwardGroup<219,7,111,28,153,102,17>(
            tables, recovery, 16, offset);
        ForwardGroup<183,224,131,222,51,187,34>(
            tables, recovery, 24, offset);
    }

    return true;
}

#undef LEO2_AVX2_T32_B256_ENTRY
#endif

}} // namespace leopard::backend
