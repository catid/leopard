/*
    Copyright (c) 2017 Christopher A. Taylor.  All rights reserved.
    See Leopard2Backend.h for the full BSD license notice.
*/

#include "Leopard2Backend.h"
#include "LeopardCommon.h"

#include <immintrin.h>

namespace leopard { namespace backend {

#if (defined(LEO2_EXPERIMENT_HIGH_T32_B256_GENERATED) || \
     defined(LEO2_EXPERIMENT_HIGH_T32_B256_TWO_BLOCK)) && \
    defined(LEO2_HAVE_AVX2_BACKEND) && !defined(NO_LEO_HAS_FF8)

namespace {

/*
    A volatile data selector keeps the control and candidate instruction
    streams identical.  Diagnostic builds change this initializer only; when
    disabled, the public routing TU falls through to the mature transform.
*/
#if defined(LEO2_EXPERIMENT_HIGH_T32_B256_GENERATED)
#ifndef LEO2_DIAGNOSTIC_DISABLE_HIGH_T32_B256_GENERATED
#define LEO2_DIAGNOSTIC_DISABLE_HIGH_T32_B256_GENERATED 0
#endif
static volatile uint32_t g_t32_b256_generated_mode =
    1U + LEO2_DIAGNOSTIC_DISABLE_HIGH_T32_B256_GENERATED;
#endif

#if defined(LEO2_EXPERIMENT_HIGH_T32_B256_TWO_BLOCK)
#ifndef LEO2_DIAGNOSTIC_DISABLE_HIGH_T32_B256_TWO_BLOCK
#define LEO2_DIAGNOSTIC_DISABLE_HIGH_T32_B256_TWO_BLOCK 0
#endif
static volatile uint32_t g_t32_b256_two_block_mode =
    1U + LEO2_DIAGNOSTIC_DISABLE_HIGH_T32_B256_TWO_BLOCK;
#endif

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
    __m256i product =
        Product(source, Broadcast(table), Broadcast(table + 16));
#if defined(__GNUC__) || defined(__clang__)
    /*
        Confine each fixed product's broadcast pair to this butterfly.  GCC
        otherwise hoists several pairs across an eight-row group, exhausting
        the sixteen AVX2 registers and spilling live codeword vectors.
    */
    /*
        Keep the product live across a compiler boundary without imposing a
        full memory-clobber scheduling barrier on every butterfly.  The read-
        write vector operand still prevents GCC from carrying the table pair
        as additional live vectors past this point.
    */
    __asm__ __volatile__("" : "+x"(product));
#endif
    destination = _mm256_xor_si256(destination, product);
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

#if defined(LEO2_EXPERIMENT_HIGH_T32_B256_TWO_BLOCK)

static LEO_FORCE_INLINE void PreparedMulAdd(
    __m256i& destination,
    __m256i source,
    __m256i low_table,
    __m256i high_table)
{
    __m256i product = Product(source, low_table, high_table);
#if defined(__GNUC__) || defined(__clang__)
    __asm__ __volatile__("" : "+x"(product));
#endif
    destination = _mm256_xor_si256(destination, product);
}

static LEO_FORCE_INLINE __m256i LoadPackedRow(
    const void* base,
    unsigned row,
    unsigned offset)
{
    const uint8_t* const bytes = static_cast<const uint8_t*>(base) +
        static_cast<size_t>(row) * 256U + offset;
    return _mm256_loadu_si256(reinterpret_cast<const __m256i*>(bytes));
}

static LEO_FORCE_INLINE void StorePackedRow(
    void* base,
    unsigned row,
    unsigned offset,
    __m256i value)
{
    uint8_t* const bytes = static_cast<uint8_t*>(base) +
        static_cast<size_t>(row) * 256U + offset;
    _mm256_storeu_si256(reinterpret_cast<__m256i*>(bytes), value);
}

static LEO_FORCE_INLINE void CandidateMemoryBoundary()
{
#if defined(__GNUC__) || defined(__clang__)
    __asm__ __volatile__("" ::: "memory");
#endif
}

template<unsigned Log>
static LEO_FORCE_INLINE void IFFT2PairPrepared(
    const uint8_t* tables,
    __m256i& x0, __m256i& y0,
    __m256i& x1, __m256i& y1)
{
    y0 = _mm256_xor_si256(y0, x0);
    y1 = _mm256_xor_si256(y1, x1);
    if (Log != 255)
    {
        const uint8_t* const table = tables + Log * 32U;
        const __m256i low = Broadcast(table);
        const __m256i high = Broadcast(table + 16);
        PreparedMulAdd(x0, y0, low, high);
        PreparedMulAdd(x1, y1, low, high);
    }
}

template<unsigned Log>
static LEO_FORCE_INLINE void IFFT2QuadPrepared(
    const uint8_t* tables,
    __m256i& x0, __m256i& y0,
    __m256i& x1, __m256i& y1,
    __m256i& x2, __m256i& y2,
    __m256i& x3, __m256i& y3)
{
    y0 = _mm256_xor_si256(y0, x0);
    y1 = _mm256_xor_si256(y1, x1);
    y2 = _mm256_xor_si256(y2, x2);
    y3 = _mm256_xor_si256(y3, x3);
    if (Log != 255)
    {
        const uint8_t* const table = tables + Log * 32U;
        const __m256i low = Broadcast(table);
        const __m256i high = Broadcast(table + 16);
        PreparedMulAdd(x0, y0, low, high);
        PreparedMulAdd(x1, y1, low, high);
        PreparedMulAdd(x2, y2, low, high);
        PreparedMulAdd(x3, y3, low, high);
    }
}

template<unsigned Log>
static LEO_FORCE_INLINE void FFT2PairPrepared(
    const uint8_t* tables,
    __m256i& x0, __m256i& y0,
    __m256i& x1, __m256i& y1)
{
    if (Log != 255)
    {
        const uint8_t* const table = tables + Log * 32U;
        const __m256i low = Broadcast(table);
        const __m256i high = Broadcast(table + 16);
        PreparedMulAdd(x0, y0, low, high);
        PreparedMulAdd(x1, y1, low, high);
    }
    y0 = _mm256_xor_si256(y0, x0);
    y1 = _mm256_xor_si256(y1, x1);
}

template<unsigned Log>
static LEO_FORCE_INLINE void FFT2QuadPrepared(
    const uint8_t* tables,
    __m256i& x0, __m256i& y0,
    __m256i& x1, __m256i& y1,
    __m256i& x2, __m256i& y2,
    __m256i& x3, __m256i& y3)
{
    if (Log != 255)
    {
        const uint8_t* const table = tables + Log * 32U;
        const __m256i low = Broadcast(table);
        const __m256i high = Broadcast(table + 16);
        PreparedMulAdd(x0, y0, low, high);
        PreparedMulAdd(x1, y1, low, high);
        PreparedMulAdd(x2, y2, low, high);
        PreparedMulAdd(x3, y3, low, high);
    }
    y0 = _mm256_xor_si256(y0, x0);
    y1 = _mm256_xor_si256(y1, x1);
    y2 = _mm256_xor_si256(y2, x2);
    y3 = _mm256_xor_si256(y3, x3);
}

template<unsigned Log01, unsigned Log23, unsigned Log02>
static LEO_FORCE_INLINE void IFFT4Prepared(
    const uint8_t* tables,
    __m256i& x0, __m256i& x1, __m256i& x2, __m256i& x3)
{
    IFFT2<Log01>(tables, x0, x1);
    IFFT2<Log23>(tables, x2, x3);
    IFFT2PairPrepared<Log02>(tables, x0, x2, x1, x3);
}

template<unsigned Log01, unsigned Log23, unsigned Log02>
static LEO_FORCE_INLINE void FFT4Prepared(
    const uint8_t* tables,
    __m256i& x0, __m256i& x1, __m256i& x2, __m256i& x3)
{
    FFT2PairPrepared<Log02>(tables, x0, x2, x1, x3);
    FFT2<Log01>(tables, x0, x1);
    FFT2<Log23>(tables, x2, x3);
}

template<unsigned Log01, unsigned Log23, unsigned Log45,
    unsigned Log67, unsigned Log02, unsigned Log46, unsigned Log04>
static LEO_FORCE_INLINE void InverseGroupPrepared(
    const uint8_t* tables,
    const void* data,
    void* work,
    unsigned base,
    unsigned offset)
{
    __m256i x0 = LoadPackedRow(data, base, offset);
    __m256i x1 = LoadPackedRow(data, base + 1, offset);
    __m256i x2 = LoadPackedRow(data, base + 2, offset);
    __m256i x3 = LoadPackedRow(data, base + 3, offset);
    __m256i x4 = LoadPackedRow(data, base + 4, offset);
    __m256i x5 = LoadPackedRow(data, base + 5, offset);
    __m256i x6 = LoadPackedRow(data, base + 6, offset);
    __m256i x7 = LoadPackedRow(data, base + 7, offset);
    IFFT2<Log01>(tables, x0, x1);
    IFFT2<Log23>(tables, x2, x3);
    IFFT2<Log45>(tables, x4, x5);
    IFFT2<Log67>(tables, x6, x7);
    IFFT2PairPrepared<Log02>(tables, x0, x2, x1, x3);
    IFFT2PairPrepared<Log46>(tables, x4, x6, x5, x7);
    IFFT2QuadPrepared<Log04>(
        tables, x0, x4, x1, x5, x2, x6, x3, x7);
    StorePackedRow(work, base, offset, x0);
    StorePackedRow(work, base + 1, offset, x1);
    StorePackedRow(work, base + 2, offset, x2);
    StorePackedRow(work, base + 3, offset, x3);
    StorePackedRow(work, base + 4, offset, x4);
    StorePackedRow(work, base + 5, offset, x5);
    StorePackedRow(work, base + 6, offset, x6);
    StorePackedRow(work, base + 7, offset, x7);
    CandidateMemoryBoundary();
}

template<unsigned Log01, unsigned Log23, unsigned Log45,
    unsigned Log67, unsigned Log02, unsigned Log46, unsigned Log04>
static LEO_FORCE_INLINE void ForwardGroupPrepared(
    const uint8_t* tables,
    void* work,
    unsigned base,
    unsigned offset)
{
    __m256i x0 = LoadPackedRow(work, base, offset);
    __m256i x1 = LoadPackedRow(work, base + 1, offset);
    __m256i x2 = LoadPackedRow(work, base + 2, offset);
    __m256i x3 = LoadPackedRow(work, base + 3, offset);
    __m256i x4 = LoadPackedRow(work, base + 4, offset);
    __m256i x5 = LoadPackedRow(work, base + 5, offset);
    __m256i x6 = LoadPackedRow(work, base + 6, offset);
    __m256i x7 = LoadPackedRow(work, base + 7, offset);
    FFT2QuadPrepared<Log04>(
        tables, x0, x4, x1, x5, x2, x6, x3, x7);
    FFT2PairPrepared<Log02>(tables, x0, x2, x1, x3);
    FFT2PairPrepared<Log46>(tables, x4, x6, x5, x7);
    FFT2<Log01>(tables, x0, x1);
    FFT2<Log23>(tables, x2, x3);
    FFT2<Log45>(tables, x4, x5);
    FFT2<Log67>(tables, x6, x7);
    StorePackedRow(work, base, offset, x0);
    StorePackedRow(work, base + 1, offset, x1);
    StorePackedRow(work, base + 2, offset, x2);
    StorePackedRow(work, base + 3, offset, x3);
    StorePackedRow(work, base + 4, offset, x4);
    StorePackedRow(work, base + 5, offset, x5);
    StorePackedRow(work, base + 6, offset, x6);
    StorePackedRow(work, base + 7, offset, x7);
    CandidateMemoryBoundary();
}

static LEO_FORCE_INLINE void InverseBlockShift32(
    const uint8_t* tables,
    const void* data,
    void* work,
    unsigned offset)
{
    InverseGroupPrepared<196,76,54,99,219,7,153>(
        tables, data, work, 0, offset);
    InverseGroupPrepared<19,49,67,52,111,28,102>(
        tables, data, work, 8, offset);
    InverseGroupPrepared<137,226,26,108,183,224,51>(
        tables, data, work, 16, offset);
    InverseGroupPrepared<27,134,38,139,131,222,187>(
        tables, data, work, 24, offset);
}

static LEO_FORCE_INLINE void InverseBlockShift64(
    const uint8_t* tables,
    const void* data,
    void* work,
    unsigned offset)
{
    InverseGroupPrepared<241,254,31,239,196,76,219>(
        tables, data, work, 0, offset);
    InverseGroupPrepared<73,200,148,140,54,99,7>(
        tables, data, work, 8, offset);
    InverseGroupPrepared<199,251,180,6,19,49,111>(
        tables, data, work, 16, offset);
    InverseGroupPrepared<35,37,2,59,67,52,28>(
        tables, data, work, 24, offset);
}

template<unsigned Log01, unsigned Log23, unsigned Log02>
static LEO_FORCE_INLINE void OuterInverse(
    const uint8_t* tables,
    void* work,
    unsigned column,
    unsigned offset)
{
    __m256i x0 = LoadPackedRow(work, column, offset);
    __m256i x1 = LoadPackedRow(work, column + 8, offset);
    __m256i x2 = LoadPackedRow(work, column + 16, offset);
    __m256i x3 = LoadPackedRow(work, column + 24, offset);
    IFFT4Prepared<Log01, Log23, Log02>(tables, x0, x1, x2, x3);
    StorePackedRow(work, column, offset, x0);
    StorePackedRow(work, column + 8, offset, x1);
    StorePackedRow(work, column + 16, offset, x2);
    StorePackedRow(work, column + 24, offset, x3);
    CandidateMemoryBoundary();
}

static LEO_FORCE_INLINE void XorStore(
    void* destination_base,
    unsigned row,
    unsigned offset,
    __m256i value)
{
    StorePackedRow(destination_base, row, offset,
        _mm256_xor_si256(
            LoadPackedRow(destination_base, row, offset), value));
}

template<unsigned Log01, unsigned Log23, unsigned Log02>
static LEO_FORCE_INLINE void OuterInverseXor(
    const uint8_t* tables,
    void* temporary,
    void* accumulator,
    unsigned column,
    unsigned offset)
{
    __m256i x0 = LoadPackedRow(temporary, column, offset);
    __m256i x1 = LoadPackedRow(temporary, column + 8, offset);
    __m256i x2 = LoadPackedRow(temporary, column + 16, offset);
    __m256i x3 = LoadPackedRow(temporary, column + 24, offset);
    IFFT4Prepared<Log01, Log23, Log02>(tables, x0, x1, x2, x3);
    XorStore(accumulator, column, offset, x0);
    XorStore(accumulator, column + 8, offset, x1);
    XorStore(accumulator, column + 16, offset, x2);
    XorStore(accumulator, column + 24, offset, x3);
    CandidateMemoryBoundary();
}

static LEO_FORCE_INLINE void OuterForward(
    const uint8_t* tables,
    void* work,
    unsigned column,
    unsigned offset)
{
    __m256i x0 = LoadPackedRow(work, column, offset);
    __m256i x1 = LoadPackedRow(work, column + 8, offset);
    __m256i x2 = LoadPackedRow(work, column + 16, offset);
    __m256i x3 = LoadPackedRow(work, column + 24, offset);
    FFT4Prepared<255,85,255>(tables, x0, x1, x2, x3);
    StorePackedRow(work, column, offset, x0);
    StorePackedRow(work, column + 8, offset, x1);
    StorePackedRow(work, column + 16, offset, x2);
    StorePackedRow(work, column + 24, offset, x3);
    CandidateMemoryBoundary();
}

static LEO_FORCE_INLINE void ForwardBlock(
    const uint8_t* tables,
    void* work,
    unsigned offset)
{
    ForwardGroupPrepared<255,85,17,34,255,85,255>(
        tables, work, 0, offset);
    ForwardGroupPrepared<153,102,51,187,17,34,85>(
        tables, work, 8, offset);
    ForwardGroupPrepared<219,7,111,28,153,102,17>(
        tables, work, 16, offset);
    ForwardGroupPrepared<183,224,131,222,51,187,34>(
        tables, work, 24, offset);
}

#endif

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

#if defined(LEO2_EXPERIMENT_HIGH_T32_B256_GENERATED)
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
#endif

#if defined(LEO2_EXPERIMENT_HIGH_T32_B256_TWO_BLOCK)
LEO2_AVX2_T32_B256_ENTRY bool TryAVX2FF8HighEncodeT32B256TwoBlockPacked(
    const void* data_base,
    void* recovery_base,
    void* temporary_base)
{
    if (g_t32_b256_two_block_mode != 1U)
        return false;

    LEO_DEBUG_ASSERT(data_base != NULL && recovery_base != NULL &&
        temporary_base != NULL);
    const uint8_t* const tables = Tables();
    LEO_DEBUG_ASSERT(tables != NULL);
    if (!tables)
        return false;

    const uint8_t* const data = static_cast<const uint8_t*>(data_base);
    for (unsigned offset = 0; offset < 256; offset += 32)
    {
        /* First message block initializes the coefficient accumulator. */
        InverseBlockShift32(tables, data, recovery_base, offset);
        for (unsigned column = 0; column < 8; ++column)
            OuterInverse<17,34,85>(tables, recovery_base, column, offset);

        /* Later blocks fold their final inverse layers into that accumulator. */
        InverseBlockShift64(
            tables, data + 32U * 256U, temporary_base, offset);
        for (unsigned column = 0; column < 8; ++column)
        {
            OuterInverseXor<153,102,17>(
                tables, temporary_base, recovery_base, column, offset);
        }

        for (unsigned column = 0; column < 8; ++column)
            OuterForward(tables, recovery_base, column, offset);
        ForwardBlock(tables, recovery_base, offset);
    }

    return true;
}
#endif

#undef LEO2_AVX2_T32_B256_ENTRY
#endif

}} // namespace leopard::backend
