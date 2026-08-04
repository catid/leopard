/*
    Copyright (c) 2017 Christopher A. Taylor.  All rights reserved.
    See Leopard2Backend.h for the full BSD license notice.
*/

#include "Leopard2Backend.h"
#include "LeopardCommon.h"

#include <cstring>
#include <immintrin.h>

namespace leopard { namespace backend {

#if defined(LEO2_HAVE_AVX2_BACKEND) && !defined(NO_LEO_HAS_FF8)

namespace {

static const size_t kNibbleHalfBytes = 16U;
static const size_t kNibbleTableBytes = 2U * kNibbleHalfBytes;
static const uint32_t kMultiply2Log = 85U;
static const uint32_t kMultiply4Log = 17U;
static_assert(kNibbleTableBytes == 32U,
    "AVX2 GF8 nibble-table record contract changed");

static LEO_FORCE_INLINE __m256i Broadcast(const uint8_t* table)
{
    return _mm256_broadcastsi128_si256(_mm_loadu_si128(
        reinterpret_cast<const __m128i*>(table)));
}

static LEO_FORCE_INLINE __m256i Product(
    __m256i data,
    __m256i low_table,
    __m256i high_table,
    __m256i mask)
{
    const __m256i low = _mm256_shuffle_epi8(
        low_table, _mm256_and_si256(data, mask));
    const __m256i high = _mm256_shuffle_epi8(
        high_table,
        _mm256_and_si256(_mm256_srli_epi64(data, 4), mask));
    return _mm256_xor_si256(low, high);
}

static LEO_FORCE_INLINE uint8_t ProductByte(
    uint8_t value,
    const uint8_t* table)
{
    return static_cast<uint8_t>(
        table[value & 15U] ^
        table[kNibbleHalfBytes + (value >> 4)]);
}

template<size_t Bytes>
struct SmallVectorIO;

template<>
struct SmallVectorIO<4U>
{
    static LEO_FORCE_INLINE __m128i Load(const uint8_t* source)
    {
        uint32_t value;
        std::memcpy(&value, source, sizeof(value));
        return _mm_cvtsi32_si128(static_cast<int>(value));
    }

    static LEO_FORCE_INLINE void Store(uint8_t* destination, __m128i value)
    {
        const uint32_t word = static_cast<uint32_t>(_mm_cvtsi128_si32(value));
        std::memcpy(destination, &word, sizeof(word));
    }
};

template<>
struct SmallVectorIO<8U>
{
    static LEO_FORCE_INLINE __m128i Load(const uint8_t* source)
    {
        return _mm_loadl_epi64(reinterpret_cast<const __m128i*>(source));
    }

    static LEO_FORCE_INLINE void Store(uint8_t* destination, __m128i value)
    {
        _mm_storel_epi64(reinterpret_cast<__m128i*>(destination), value);
    }
};

template<>
struct SmallVectorIO<16U>
{
    static LEO_FORCE_INLINE __m128i Load(const uint8_t* source)
    {
        return _mm_loadu_si128(reinterpret_cast<const __m128i*>(source));
    }

    static LEO_FORCE_INLINE void Store(uint8_t* destination, __m128i value)
    {
        _mm_storeu_si128(reinterpret_cast<__m128i*>(destination), value);
    }
};

static LEO_FORCE_INLINE __m128i Product128(
    __m128i data,
    __m128i low_table,
    __m128i high_table,
    __m128i mask)
{
    const __m128i low = _mm_shuffle_epi8(
        low_table, _mm_and_si128(data, mask));
    const __m128i high = _mm_shuffle_epi8(
        high_table,
        _mm_and_si128(_mm_srli_epi64(data, 4), mask));
    return _mm_xor_si128(low, high);
}

template<size_t Bytes>
static LEO_FORCE_INLINE void EncodeSmallVector(
    const uint8_t* const* input,
    uint8_t* output0,
    uint8_t* output1,
    uint64_t offset,
    __m128i multiply2_low,
    __m128i multiply2_high,
    __m128i multiply4_low,
    __m128i multiply4_high,
    __m128i mask)
{
    const __m128i a = SmallVectorIO<Bytes>::Load(input[0] + offset);
    const __m128i b = SmallVectorIO<Bytes>::Load(input[1] + offset);
    const __m128i c = SmallVectorIO<Bytes>::Load(input[2] + offset);
    const __m128i d = SmallVectorIO<Bytes>::Load(input[3] + offset);
    const __m128i common = _mm_xor_si128(
        Product128(_mm_xor_si128(a, b),
            multiply2_low, multiply2_high, mask),
        Product128(_mm_xor_si128(c, d),
            multiply4_low, multiply4_high, mask));
    const __m128i parity0 = _mm_xor_si128(
        _mm_xor_si128(a, c), common);
    const __m128i parity1 = _mm_xor_si128(
        _mm_xor_si128(b, d), common);
    SmallVectorIO<Bytes>::Store(output0 + offset, parity0);
    SmallVectorIO<Bytes>::Store(output1 + offset, parity1);
}

static LEO_FORCE_INLINE void EncodeVector(
    const uint8_t* const* input,
    uint8_t* output0,
    uint8_t* output1,
    uint64_t offset,
    __m256i multiply2_low,
    __m256i multiply2_high,
    __m256i multiply4_low,
    __m256i multiply4_high,
    __m256i mask)
{
    const __m256i a = _mm256_loadu_si256(
        reinterpret_cast<const __m256i*>(input[0] + offset));
    const __m256i b = _mm256_loadu_si256(
        reinterpret_cast<const __m256i*>(input[1] + offset));
    const __m256i c = _mm256_loadu_si256(
        reinterpret_cast<const __m256i*>(input[2] + offset));
    const __m256i d = _mm256_loadu_si256(
        reinterpret_cast<const __m256i*>(input[3] + offset));
    const __m256i common = _mm256_xor_si256(
        Product(_mm256_xor_si256(a, b),
            multiply2_low, multiply2_high, mask),
        Product(_mm256_xor_si256(c, d),
            multiply4_low, multiply4_high, mask));
    const __m256i parity0 = _mm256_xor_si256(
        _mm256_xor_si256(a, c), common);
    const __m256i parity1 = _mm256_xor_si256(
        _mm256_xor_si256(b, d), common);
    _mm256_storeu_si256(
        reinterpret_cast<__m256i*>(output0 + offset), parity0);
    _mm256_storeu_si256(
        reinterpret_cast<__m256i*>(output1 + offset), parity1);
}

#if defined(_MSC_VER)
#define LEO2_AVX2_T2_K4_ENTRY __declspec(noinline)
#elif defined(__GNUC__) && !defined(__clang__) && defined(__ELF__)
#define LEO2_AVX2_T2_K4_ENTRY \
    __attribute__((noinline, noipa, section(".text.leo2_t2_k4"), aligned(64)))
#elif defined(__clang__) && defined(__ELF__)
#define LEO2_AVX2_T2_K4_ENTRY \
    __attribute__((noinline, section(".text.leo2_t2_k4"), aligned(64)))
#elif defined(__GNUC__) || defined(__clang__)
#define LEO2_AVX2_T2_K4_ENTRY __attribute__((noinline, aligned(64)))
#else
#define LEO2_AVX2_T2_K4_ENTRY
#endif

} // namespace

LEO2_AVX2_T2_K4_ENTRY void AVX2FF8HighEncodeT2K4Packed(
    const void* const* data,
    void* const* recovery,
    uint64_t byte_count)
{
    LEO_DEBUG_ASSERT(data != NULL && recovery != NULL && byte_count != 0);
    const uint8_t* const tables = static_cast<const uint8_t*>(
        GetAVX2FF8Tables());
    LEO_DEBUG_ASSERT(tables != NULL);

    /* The independent systematic generator rows are

           p0 = 3*a + 2*b + 5*c + 4*d
           p1 = 2*a + 3*b + 4*c + 5*d.

       In characteristic two this is one shared pair of fixed products:

           q = 2*(a+b) + 4*(c+d)
           p0 = a+c+q, p1 = b+d+q.

       Leopard's nibble-table indices are logarithms, with log(2)=85 and
       log(4)=17 in the legacy Cantor representation. */
    const uint8_t* const multiply2 =
        tables + kMultiply2Log * kNibbleTableBytes;
    const uint8_t* const multiply4 =
        tables + kMultiply4Log * kNibbleTableBytes;
    const uint8_t* input[4] = {
        static_cast<const uint8_t*>(data[0]),
        static_cast<const uint8_t*>(data[1]),
        static_cast<const uint8_t*>(data[2]),
        static_cast<const uint8_t*>(data[3])
    };
    uint8_t* const output0 = static_cast<uint8_t*>(recovery[0]);
    uint8_t* const output1 = static_cast<uint8_t*>(recovery[1]);

    uint64_t offset = 0;
    if (byte_count >= 32U)
    {
        /* These four fixed tables are invariant across every byte chunk.
           Keeping them in YMM registers removes four broadcasts per loop and
           also supplies the narrower overlapping final-vector paths. */
        const __m256i mask = _mm256_set1_epi8(15);
        const __m256i multiply2_low = Broadcast(multiply2);
        const __m256i multiply2_high =
            Broadcast(multiply2 + kNibbleHalfBytes);
        const __m256i multiply4_low = Broadcast(multiply4);
        const __m256i multiply4_high =
            Broadcast(multiply4 + kNibbleHalfBytes);
        do
        {
            EncodeVector(input, output0, output1, offset,
                multiply2_low, multiply2_high,
                multiply4_low, multiply4_high, mask);
            offset += 32U;
        }
        while (byte_count - offset >= 32U);

        const uint64_t tail_bytes = byte_count - offset;
        if (tail_bytes == 0)
            return;

        /* Every byte is independent, so a final vector may safely overlap
           the preceding vector.  For a 16..31-byte remainder this keeps the
           load/store interval wholly inside [0, byte_count) while avoiding
           the substantially slower scalar nibble-table loop. */
        if (tail_bytes >= 16U)
        {
            EncodeVector(input, output0, output1, byte_count - 32U,
                multiply2_low, multiply2_high,
                multiply4_low, multiply4_high, mask);
            return;
        }

        const __m128i mask128 = _mm256_castsi256_si128(mask);
        const __m128i multiply2_low128 =
            _mm256_castsi256_si128(multiply2_low);
        const __m128i multiply2_high128 =
            _mm256_castsi256_si128(multiply2_high);
        const __m128i multiply4_low128 =
            _mm256_castsi256_si128(multiply4_low);
        const __m128i multiply4_high128 =
            _mm256_castsi256_si128(multiply4_high);
        if (tail_bytes >= 8U)
        {
            EncodeSmallVector<16U>(input, output0, output1, byte_count - 16U,
                multiply2_low128, multiply2_high128,
                multiply4_low128, multiply4_high128, mask128);
            return;
        }
        if (tail_bytes >= 4U)
        {
            EncodeSmallVector<8U>(input, output0, output1, byte_count - 8U,
                multiply2_low128, multiply2_high128,
                multiply4_low128, multiply4_high128, mask128);
            return;
        }
        EncodeSmallVector<4U>(input, output0, output1, byte_count - 4U,
            multiply2_low128, multiply2_high128,
            multiply4_low128, multiply4_high128, mask128);
        return;
    }
    else if (byte_count >= 16U)
    {
        const __m128i mask = _mm_set1_epi8(15);
        const __m128i multiply2_low = _mm_loadu_si128(
            reinterpret_cast<const __m128i*>(multiply2));
        const __m128i multiply2_high = _mm_loadu_si128(
            reinterpret_cast<const __m128i*>(
                multiply2 + kNibbleHalfBytes));
        const __m128i multiply4_low = _mm_loadu_si128(
            reinterpret_cast<const __m128i*>(multiply4));
        const __m128i multiply4_high = _mm_loadu_si128(
            reinterpret_cast<const __m128i*>(
                multiply4 + kNibbleHalfBytes));
        EncodeSmallVector<16U>(input, output0, output1, 0,
            multiply2_low, multiply2_high,
            multiply4_low, multiply4_high, mask);
        if (byte_count != 16U)
        {
            EncodeSmallVector<16U>(input, output0, output1, byte_count - 16U,
                multiply2_low, multiply2_high,
                multiply4_low, multiply4_high, mask);
        }
        return;
    }
    else if (byte_count >= 8U)
    {
        const __m128i mask = _mm_set1_epi8(15);
        const __m128i multiply2_low = _mm_loadu_si128(
            reinterpret_cast<const __m128i*>(multiply2));
        const __m128i multiply2_high = _mm_loadu_si128(
            reinterpret_cast<const __m128i*>(
                multiply2 + kNibbleHalfBytes));
        const __m128i multiply4_low = _mm_loadu_si128(
            reinterpret_cast<const __m128i*>(multiply4));
        const __m128i multiply4_high = _mm_loadu_si128(
            reinterpret_cast<const __m128i*>(
                multiply4 + kNibbleHalfBytes));
        EncodeSmallVector<8U>(input, output0, output1, 0,
            multiply2_low, multiply2_high,
            multiply4_low, multiply4_high, mask);
        if (byte_count != 8U)
        {
            EncodeSmallVector<8U>(input, output0, output1, byte_count - 8U,
                multiply2_low, multiply2_high,
                multiply4_low, multiply4_high, mask);
        }
        return;
    }
    for (; offset < byte_count; ++offset)
    {
        const uint8_t a = input[0][offset];
        const uint8_t b = input[1][offset];
        const uint8_t c = input[2][offset];
        const uint8_t d = input[3][offset];
        const uint8_t common = static_cast<uint8_t>(
            ProductByte(static_cast<uint8_t>(a ^ b), multiply2) ^
            ProductByte(static_cast<uint8_t>(c ^ d), multiply4));
        output0[offset] = static_cast<uint8_t>(a ^ c ^ common);
        output1[offset] = static_cast<uint8_t>(b ^ d ^ common);
    }
}

#undef LEO2_AVX2_T2_K4_ENTRY

#endif

}} // namespace leopard::backend
