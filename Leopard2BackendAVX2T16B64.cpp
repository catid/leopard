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

#endif

}} // namespace leopard::backend
