/*
    Copyright (c) 2017 Christopher A. Taylor.  All rights reserved.
    See Leopard2Backend.h for the full BSD license notice.
*/

#include "Leopard2Backend.h"
#include "LeopardCommon.h"

#include <immintrin.h>

namespace leopard { namespace backend {

#if defined(LEO2_EXPERIMENT_HIGH_T16_B64_GENERATED) && \
    defined(LEO2_HAVE_AVX2_BACKEND) && !defined(NO_LEO_HAS_FF8)

namespace {

/*
    Generated exact K=R=T=16, B=64 legacy-high transform.

    The ordinary transform makes four complete shard passes: contiguous
    inverse radix-four groups, an outer inverse group, the matching outer
    forward group, and contiguous forward groups.  The two outer groups touch
    the same four columns, so this kernel keeps each column in registers and
    executes both groups before storing it.  This removes one 16-shard memory
    round while preserving the legacy butterfly order and parity bytes.

    Keep this code in a trailing AVX2-only object.  Earlier in-TU experiments
    moved mature hot functions enough to regress unrelated T=32 tiny-shard
    encoding even though the generated arithmetic itself was substantially
    faster.  The ordinary AVX2 backend owns the immutable nibble tables; this
    object reads their documented byte representation.
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

static LEO_FORCE_INLINE __m256i Load(
    const void* pointer,
    unsigned offset)
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
    static_assert(Log <= 255, "T16 multiplier log is out of range");
    const uint8_t* const table = tables + Log * 32U;
    __m256i product = Product(
        source, Broadcast(table), Broadcast(table + 16));
#if defined(__GNUC__) || defined(__clang__)
    /* Bound each table pair's live range on the sixteen-register AVX2 file. */
    __asm__ __volatile__("" : "+x"(product));
#endif
    destination = _mm256_xor_si256(destination, product);
}

template<unsigned Log01, unsigned Log23, unsigned Log02>
static LEO_FORCE_INLINE void IFFT4(
    const uint8_t* tables,
    __m256i& value0,
    __m256i& value1,
    __m256i& value2,
    __m256i& value3)
{
    value1 = _mm256_xor_si256(value1, value0);
    MulAdd<Log01>(tables, value0, value1);
    value3 = _mm256_xor_si256(value3, value2);
    MulAdd<Log23>(tables, value2, value3);
    value2 = _mm256_xor_si256(value2, value0);
    value3 = _mm256_xor_si256(value3, value1);
    MulAdd<Log02>(tables, value0, value2);
    MulAdd<Log02>(tables, value1, value3);
}

template<unsigned Log01, unsigned Log23, unsigned Log02>
static LEO_FORCE_INLINE void FFT4(
    const uint8_t* tables,
    __m256i& value0,
    __m256i& value1,
    __m256i& value2,
    __m256i& value3)
{
    if (Log02 != 255)
    {
        MulAdd<Log02>(tables, value0, value2);
        MulAdd<Log02>(tables, value1, value3);
    }
    value2 = _mm256_xor_si256(value2, value0);
    value3 = _mm256_xor_si256(value3, value1);
    if (Log01 != 255)
        MulAdd<Log01>(tables, value0, value1);
    value1 = _mm256_xor_si256(value1, value0);
    if (Log23 != 255)
        MulAdd<Log23>(tables, value2, value3);
    value3 = _mm256_xor_si256(value3, value2);
}

template<unsigned Log01, unsigned Log23, unsigned Log02>
static LEO_FORCE_INLINE void InverseGroup(
    const uint8_t* tables,
    const void* const* data,
    void* const* work,
    unsigned base,
    unsigned offset)
{
    __m256i value0 = Load(data[base], offset);
    __m256i value1 = Load(data[base + 1], offset);
    __m256i value2 = Load(data[base + 2], offset);
    __m256i value3 = Load(data[base + 3], offset);
    IFFT4<Log01, Log23, Log02>(
        tables, value0, value1, value2, value3);
    Store(work[base], offset, value0);
    Store(work[base + 1], offset, value1);
    Store(work[base + 2], offset, value2);
    Store(work[base + 3], offset, value3);
}

static LEO_FORCE_INLINE void OuterGroups(
    const uint8_t* tables,
    void* const* work,
    unsigned column,
    unsigned offset)
{
    __m256i value0 = Load(work[column], offset);
    __m256i value1 = Load(work[column + 4], offset);
    __m256i value2 = Load(work[column + 8], offset);
    __m256i value3 = Load(work[column + 12], offset);
    IFFT4<17, 34, 85>(tables, value0, value1, value2, value3);
    FFT4<255, 85, 255>(tables, value0, value1, value2, value3);
    Store(work[column], offset, value0);
    Store(work[column + 4], offset, value1);
    Store(work[column + 8], offset, value2);
    Store(work[column + 12], offset, value3);
}

template<unsigned Log01, unsigned Log23, unsigned Log02>
static LEO_FORCE_INLINE void ForwardGroup(
    const uint8_t* tables,
    void* const* work,
    unsigned base,
    unsigned offset)
{
    __m256i value0 = Load(work[base], offset);
    __m256i value1 = Load(work[base + 1], offset);
    __m256i value2 = Load(work[base + 2], offset);
    __m256i value3 = Load(work[base + 3], offset);
    FFT4<Log01, Log23, Log02>(
        tables, value0, value1, value2, value3);
    Store(work[base], offset, value0);
    Store(work[base + 1], offset, value1);
    Store(work[base + 2], offset, value2);
    Store(work[base + 3], offset, value3);
}

#if defined(_MSC_VER)
#define LEO2_AVX2_T16_B64_ENTRY __declspec(noinline)
#elif defined(__GNUC__) && !defined(__clang__) && defined(__ELF__)
#define LEO2_AVX2_T16_B64_ENTRY \
    __attribute__((noinline, noipa, section(".text.leo2_t16_b64"), aligned(64)))
#elif defined(__clang__) && defined(__ELF__)
#define LEO2_AVX2_T16_B64_ENTRY \
    __attribute__((noinline, section(".text.leo2_t16_b64"), aligned(64)))
#elif defined(__GNUC__) || defined(__clang__)
#define LEO2_AVX2_T16_B64_ENTRY __attribute__((noinline, aligned(64)))
#else
#define LEO2_AVX2_T16_B64_ENTRY
#endif

static LEO2_AVX2_T16_B64_ENTRY void Encode(
    const uint8_t* tables,
    const void* const* data,
    void* const* work)
{
    InverseGroup<219, 7, 153>(tables, data, work, 0, 0);
    InverseGroup<111, 28, 102>(tables, data, work, 4, 0);
    InverseGroup<183, 224, 51>(tables, data, work, 8, 0);
    InverseGroup<131, 222, 187>(tables, data, work, 12, 0);
    InverseGroup<219, 7, 153>(tables, data, work, 0, 32);
    InverseGroup<111, 28, 102>(tables, data, work, 4, 32);
    InverseGroup<183, 224, 51>(tables, data, work, 8, 32);
    InverseGroup<131, 222, 187>(tables, data, work, 12, 32);

    OuterGroups(tables, work, 0, 0);
    OuterGroups(tables, work, 1, 0);
    OuterGroups(tables, work, 2, 0);
    OuterGroups(tables, work, 3, 0);
    OuterGroups(tables, work, 0, 32);
    OuterGroups(tables, work, 1, 32);
    OuterGroups(tables, work, 2, 32);
    OuterGroups(tables, work, 3, 32);

    ForwardGroup<255, 85, 255>(tables, work, 0, 0);
    ForwardGroup<17, 34, 85>(tables, work, 4, 0);
    ForwardGroup<153, 102, 17>(tables, work, 8, 0);
    ForwardGroup<51, 187, 34>(tables, work, 12, 0);
    ForwardGroup<255, 85, 255>(tables, work, 0, 32);
    ForwardGroup<17, 34, 85>(tables, work, 4, 32);
    ForwardGroup<153, 102, 17>(tables, work, 8, 32);
    ForwardGroup<51, 187, 34>(tables, work, 12, 32);
}

struct PreparedProduct
{
    __m256i low;
    __m256i high;
};

template<unsigned Log>
static LEO_FORCE_INLINE PreparedProduct PrepareProduct(
    const uint8_t* tables)
{
    static_assert(Log <= 255, "T16 prepared multiplier log is out of range");
    PreparedProduct product;
    if (Log == 255)
    {
        product.low = _mm256_setzero_si256();
        product.high = _mm256_setzero_si256();
    }
    else
    {
        const uint8_t* const table = tables + Log * 32U;
        product.low = Broadcast(table);
        product.high = Broadcast(table + 16);
    }
    return product;
}

static LEO_FORCE_INLINE void MulAddPrepared(
    __m256i& destination,
    __m256i source,
    const PreparedProduct& product)
{
    destination = _mm256_xor_si256(destination,
        Product(source, product.low, product.high));
}

static LEO_FORCE_INLINE void IFFT4Prepared(
    __m256i& value0,
    __m256i& value1,
    __m256i& value2,
    __m256i& value3,
    const PreparedProduct& product01,
    const PreparedProduct& product23,
    const PreparedProduct& product02)
{
    value1 = _mm256_xor_si256(value1, value0);
    MulAddPrepared(value0, value1, product01);
    value3 = _mm256_xor_si256(value3, value2);
    MulAddPrepared(value2, value3, product23);
    value2 = _mm256_xor_si256(value2, value0);
    value3 = _mm256_xor_si256(value3, value1);
    MulAddPrepared(value0, value2, product02);
    MulAddPrepared(value1, value3, product02);
}

template<bool Use01, bool Use23, bool Use02>
static LEO_FORCE_INLINE void FFT4Prepared(
    __m256i& value0,
    __m256i& value1,
    __m256i& value2,
    __m256i& value3,
    const PreparedProduct& product01,
    const PreparedProduct& product23,
    const PreparedProduct& product02)
{
    if (Use02)
    {
        MulAddPrepared(value0, value2, product02);
        MulAddPrepared(value1, value3, product02);
    }
    value2 = _mm256_xor_si256(value2, value0);
    value3 = _mm256_xor_si256(value3, value1);
    if (Use01)
        MulAddPrepared(value0, value1, product01);
    value1 = _mm256_xor_si256(value1, value0);
    if (Use23)
        MulAddPrepared(value2, value3, product23);
    value3 = _mm256_xor_si256(value3, value2);
}

template<unsigned Log01, unsigned Log23, unsigned Log02>
static LEO_FORCE_INLINE void InverseGroupPrepared(
    const uint8_t* tables,
    const void* const* data,
    void* const* work,
    unsigned base,
    unsigned byte_count)
{
    const PreparedProduct product01 = PrepareProduct<Log01>(tables);
    const PreparedProduct product23 = PrepareProduct<Log23>(tables);
    const PreparedProduct product02 = PrepareProduct<Log02>(tables);
    for (unsigned offset = 0; offset < byte_count; offset += 32)
    {
        __m256i value0 = Load(data[base], offset);
        __m256i value1 = Load(data[base + 1], offset);
        __m256i value2 = Load(data[base + 2], offset);
        __m256i value3 = Load(data[base + 3], offset);
        IFFT4Prepared(value0, value1, value2, value3,
            product01, product23, product02);
        Store(work[base], offset, value0);
        Store(work[base + 1], offset, value1);
        Store(work[base + 2], offset, value2);
        Store(work[base + 3], offset, value3);
    }
}

static LEO_FORCE_INLINE void OuterGroupsPrepared(
    const uint8_t* tables,
    void* const* work,
    unsigned byte_count)
{
    const PreparedProduct product17 = PrepareProduct<17>(tables);
    const PreparedProduct product34 = PrepareProduct<34>(tables);
    const PreparedProduct product85 = PrepareProduct<85>(tables);
    for (unsigned column = 0; column < 4; ++column)
    {
        for (unsigned offset = 0; offset < byte_count; offset += 32)
        {
            __m256i value0 = Load(work[column], offset);
            __m256i value1 = Load(work[column + 4], offset);
            __m256i value2 = Load(work[column + 8], offset);
            __m256i value3 = Load(work[column + 12], offset);
            IFFT4Prepared(value0, value1, value2, value3,
                product17, product34, product85);
            FFT4Prepared<false, true, false>(
                value0, value1, value2, value3,
                product85, product85, product85);
            Store(work[column], offset, value0);
            Store(work[column + 4], offset, value1);
            Store(work[column + 8], offset, value2);
            Store(work[column + 12], offset, value3);
        }
    }
}

template<unsigned Log01, unsigned Log23, unsigned Log02>
static LEO_FORCE_INLINE void ForwardGroupPrepared(
    const uint8_t* tables,
    void* const* work,
    unsigned base,
    unsigned byte_count)
{
    const PreparedProduct product01 = PrepareProduct<Log01>(tables);
    const PreparedProduct product23 = PrepareProduct<Log23>(tables);
    const PreparedProduct product02 = PrepareProduct<Log02>(tables);
    for (unsigned offset = 0; offset < byte_count; offset += 32)
    {
        __m256i value0 = Load(work[base], offset);
        __m256i value1 = Load(work[base + 1], offset);
        __m256i value2 = Load(work[base + 2], offset);
        __m256i value3 = Load(work[base + 3], offset);
        FFT4Prepared<Log01 != 255, Log23 != 255, Log02 != 255>(
            value0, value1, value2, value3,
            product01, product23, product02);
        Store(work[base], offset, value0);
        Store(work[base + 1], offset, value1);
        Store(work[base + 2], offset, value2);
        Store(work[base + 3], offset, value3);
    }
}

static LEO2_AVX2_T16_B64_ENTRY void EncodePrepared(
    const uint8_t* tables,
    const void* const* data,
    void* const* work,
    unsigned byte_count)
{
    InverseGroupPrepared<219, 7, 153>(
        tables, data, work, 0, byte_count);
    InverseGroupPrepared<111, 28, 102>(
        tables, data, work, 4, byte_count);
    InverseGroupPrepared<183, 224, 51>(
        tables, data, work, 8, byte_count);
    InverseGroupPrepared<131, 222, 187>(
        tables, data, work, 12, byte_count);
    OuterGroupsPrepared(tables, work, byte_count);
    ForwardGroupPrepared<255, 85, 255>(
        tables, work, 0, byte_count);
    ForwardGroupPrepared<17, 34, 85>(
        tables, work, 4, byte_count);
    ForwardGroupPrepared<153, 102, 17>(
        tables, work, 8, byte_count);
    ForwardGroupPrepared<51, 187, 34>(
        tables, work, 12, byte_count);
}

template<unsigned OriginalCount>
static LEO_FORCE_INLINE void InverseTailPrepared(
    const uint8_t* tables,
    const void* const* data,
    void* const* work,
    unsigned byte_count)
{
    static_assert(OriginalCount >= 9 && OriginalCount <= 12,
        "T16 inverse-tail count is outside the partial group");
    const PreparedProduct product01 = PrepareProduct<183>(tables);
    const PreparedProduct product23 = PrepareProduct<224>(tables);
    const PreparedProduct product02 = PrepareProduct<51>(tables);
    for (unsigned offset = 0; offset < byte_count; offset += 32)
    {
        __m256i value0 = Load(data[8], offset);
        __m256i value1 = OriginalCount >= 10
            ? Load(data[9], offset) : _mm256_setzero_si256();
        __m256i value2 = OriginalCount >= 11
            ? Load(data[10], offset) : _mm256_setzero_si256();
        __m256i value3 = OriginalCount >= 12
            ? Load(data[11], offset) : _mm256_setzero_si256();
        value1 = _mm256_xor_si256(value1, value0);
        MulAddPrepared(value0, value1, product01);
        if (OriginalCount >= 11)
        {
            value3 = _mm256_xor_si256(value3, value2);
            MulAddPrepared(value2, value3, product23);
        }
        value2 = _mm256_xor_si256(value2, value0);
        value3 = _mm256_xor_si256(value3, value1);
        MulAddPrepared(value0, value2, product02);
        MulAddPrepared(value1, value3, product02);
        Store(work[8], offset, value0);
        Store(work[9], offset, value1);
        Store(work[10], offset, value2);
        Store(work[11], offset, value3);
    }
}

static LEO_FORCE_INLINE void OuterGroupsPreparedLastZero(
    const uint8_t* tables,
    void* const* work,
    unsigned byte_count)
{
    const PreparedProduct product17 = PrepareProduct<17>(tables);
    const PreparedProduct product34 = PrepareProduct<34>(tables);
    const PreparedProduct product85 = PrepareProduct<85>(tables);
    for (unsigned column = 0; column < 4; ++column)
    {
        for (unsigned offset = 0; offset < byte_count; offset += 32)
        {
            __m256i value0 = Load(work[column], offset);
            __m256i value1 = Load(work[column + 4], offset);
            __m256i value2 = Load(work[column + 8], offset);
            __m256i value3 = _mm256_setzero_si256();
            IFFT4Prepared(value0, value1, value2, value3,
                product17, product34, product85);

            /* FFT4<255,85,255>, omitting only the unused fourth output. */
            value2 = _mm256_xor_si256(value2, value0);
            value3 = _mm256_xor_si256(value3, value1);
            value1 = _mm256_xor_si256(value1, value0);
            MulAddPrepared(value2, value3, product85);

            Store(work[column], offset, value0);
            Store(work[column + 4], offset, value1);
            Store(work[column + 8], offset, value2);
        }
    }
}

template<unsigned OriginalCount>
static LEO_FORCE_INLINE void ForwardTailPrepared(
    const uint8_t* tables,
    void* const* work,
    unsigned byte_count)
{
    static_assert(OriginalCount >= 9 && OriginalCount <= 12,
        "T16 forward-tail count is outside the partial group");
    const PreparedProduct product01 = PrepareProduct<153>(tables);
    const PreparedProduct product23 = PrepareProduct<102>(tables);
    const PreparedProduct product02 = PrepareProduct<17>(tables);
    for (unsigned offset = 0; offset < byte_count; offset += 32)
    {
        __m256i value0 = Load(work[8], offset);
        __m256i value1 = Load(work[9], offset);
        __m256i value2 = Load(work[10], offset);
        __m256i value3 = Load(work[11], offset);
        MulAddPrepared(value0, value2, product02);
        MulAddPrepared(value1, value3, product02);
        if (OriginalCount >= 11)
        {
            value2 = _mm256_xor_si256(value2, value0);
            value3 = _mm256_xor_si256(value3, value1);
        }
        MulAddPrepared(value0, value1, product01);
        value1 = _mm256_xor_si256(value1, value0);
        if (OriginalCount >= 11)
            MulAddPrepared(value2, value3, product23);
        if (OriginalCount >= 12)
            value3 = _mm256_xor_si256(value3, value2);
        Store(work[8], offset, value0);
        if (OriginalCount >= 10)
            Store(work[9], offset, value1);
        if (OriginalCount >= 11)
            Store(work[10], offset, value2);
        if (OriginalCount >= 12)
            Store(work[11], offset, value3);
    }
}

template<unsigned OriginalCount>
static LEO_FORCE_INLINE void InverseFinalPrepared(
    const uint8_t* tables,
    const void* const* data,
    void* const* work,
    unsigned byte_count)
{
    static_assert(OriginalCount >= 13 && OriginalCount <= 15,
        "T16 final inverse count is outside the partial group");
    const PreparedProduct product01 = PrepareProduct<131>(tables);
    const PreparedProduct product23 = PrepareProduct<222>(tables);
    const PreparedProduct product02 = PrepareProduct<187>(tables);
    for (unsigned offset = 0; offset < byte_count; offset += 32)
    {
        __m256i value0 = Load(data[12], offset);
        __m256i value1 = OriginalCount >= 14
            ? Load(data[13], offset) : _mm256_setzero_si256();
        __m256i value2 = OriginalCount >= 15
            ? Load(data[14], offset) : _mm256_setzero_si256();
        __m256i value3 = _mm256_setzero_si256();
        value1 = _mm256_xor_si256(value1, value0);
        MulAddPrepared(value0, value1, product01);
        if (OriginalCount >= 15)
        {
            value3 = value2;
            MulAddPrepared(value2, value3, product23);
        }
        value2 = _mm256_xor_si256(value2, value0);
        value3 = _mm256_xor_si256(value3, value1);
        MulAddPrepared(value0, value2, product02);
        MulAddPrepared(value1, value3, product02);
        Store(work[12], offset, value0);
        Store(work[13], offset, value1);
        Store(work[14], offset, value2);
        Store(work[15], offset, value3);
    }
}

template<unsigned OriginalCount>
static LEO_FORCE_INLINE void ForwardFinalPrepared(
    const uint8_t* tables,
    void* const* work,
    unsigned byte_count)
{
    static_assert(OriginalCount >= 13 && OriginalCount <= 15,
        "T16 final forward count is outside the partial group");
    const PreparedProduct product01 = PrepareProduct<51>(tables);
    const PreparedProduct product23 = PrepareProduct<187>(tables);
    const PreparedProduct product02 = PrepareProduct<34>(tables);
    for (unsigned offset = 0; offset < byte_count; offset += 32)
    {
        __m256i value0 = Load(work[12], offset);
        __m256i value1 = Load(work[13], offset);
        __m256i value2 = Load(work[14], offset);
        __m256i value3 = Load(work[15], offset);
        MulAddPrepared(value0, value2, product02);
        MulAddPrepared(value1, value3, product02);
        if (OriginalCount >= 15)
        {
            value2 = _mm256_xor_si256(value2, value0);
            value3 = _mm256_xor_si256(value3, value1);
        }
        MulAddPrepared(value0, value1, product01);
        value1 = _mm256_xor_si256(value1, value0);
        if (OriginalCount >= 15)
            MulAddPrepared(value2, value3, product23);
        Store(work[12], offset, value0);
        if (OriginalCount >= 14)
            Store(work[13], offset, value1);
        if (OriginalCount >= 15)
            Store(work[14], offset, value2);
    }
}

template<unsigned OriginalCount>
static LEO2_AVX2_T16_B64_ENTRY void EncodePreparedPartial(
    const uint8_t* tables,
    const void* const* data,
    void* const* work,
    unsigned byte_count)
{
    InverseGroupPrepared<219, 7, 153>(
        tables, data, work, 0, byte_count);
    InverseGroupPrepared<111, 28, 102>(
        tables, data, work, 4, byte_count);
    InverseTailPrepared<OriginalCount>(
        tables, data, work, byte_count);
    OuterGroupsPreparedLastZero(tables, work, byte_count);
    ForwardGroupPrepared<255, 85, 255>(
        tables, work, 0, byte_count);
    ForwardGroupPrepared<17, 34, 85>(
        tables, work, 4, byte_count);
    ForwardTailPrepared<OriginalCount>(tables, work, byte_count);
}

static LEO2_AVX2_T16_B64_ENTRY void EncodePreparedFinalPartial(
    const uint8_t* tables,
    const void* const* data,
    void* const* work,
    unsigned original_count,
    unsigned byte_count)
{
    InverseGroupPrepared<219, 7, 153>(
        tables, data, work, 0, byte_count);
    InverseGroupPrepared<111, 28, 102>(
        tables, data, work, 4, byte_count);
    InverseGroupPrepared<183, 224, 51>(
        tables, data, work, 8, byte_count);
    switch (original_count)
    {
    case 13: InverseFinalPrepared<13>(tables, data, work, byte_count); break;
    case 14: InverseFinalPrepared<14>(tables, data, work, byte_count); break;
    default: InverseFinalPrepared<15>(tables, data, work, byte_count); break;
    }
    OuterGroupsPrepared(tables, work, byte_count);
    ForwardGroupPrepared<255, 85, 255>(
        tables, work, 0, byte_count);
    ForwardGroupPrepared<17, 34, 85>(
        tables, work, 4, byte_count);
    ForwardGroupPrepared<153, 102, 17>(
        tables, work, 8, byte_count);
    switch (original_count)
    {
    case 13: ForwardFinalPrepared<13>(tables, work, byte_count); break;
    case 14: ForwardFinalPrepared<14>(tables, work, byte_count); break;
    default: ForwardFinalPrepared<15>(tables, work, byte_count); break;
    }
}

#undef LEO2_AVX2_T16_B64_ENTRY

} // namespace

bool TryAVX2FF8HighEncodeT16B64(
    const void* const* data,
    void* const* recovery)
{
    const uint8_t* const tables = Tables();
    if (!tables || !data || !recovery)
        return false;
    Encode(tables, data, recovery);
    return true;
}

bool TryAVX2FF8HighEncodeT16(
    const void* const* data,
    void* const* recovery,
    uint32_t original_count,
    uint64_t byte_count)
{
    const uint8_t* const tables = Tables();
    if (!tables || !data || !recovery || byte_count == 0 ||
        byte_count > UINT32_MAX || (byte_count & 63U) != 0 ||
        original_count < 9 || original_count > 16)
        return false;
    if (byte_count == 64 && original_count == 16)
    {
        Encode(tables, data, recovery);
        return true;
    }
    const unsigned bytes = static_cast<unsigned>(byte_count);
    switch (original_count)
    {
    case 9: EncodePreparedPartial<9>(tables, data, recovery, bytes); break;
    case 10: EncodePreparedPartial<10>(tables, data, recovery, bytes); break;
    case 11: EncodePreparedPartial<11>(tables, data, recovery, bytes); break;
    case 12: EncodePreparedPartial<12>(tables, data, recovery, bytes); break;
    case 13:
    case 14:
    case 15:
        EncodePreparedFinalPartial(
            tables, data, recovery, original_count, bytes);
        break;
    default: EncodePrepared(tables, data, recovery, bytes); break;
    }
    return true;
}

#endif

}} // namespace leopard::backend
