/*
    Copyright (c) 2017 Christopher A. Taylor.  All rights reserved.
    See Leopard2Backend.h for the full BSD license notice.
*/

#include "Leopard2Backend.h"
#include "LeopardCommon.h"

#include <cstddef>
#include <cstdint>
#include <cstring>

#include <immintrin.h>

#if !defined(LEO2_HAVE_AVX512_GFNI_T16)
#error "This translation unit requires LEO2_HAVE_AVX512_GFNI_T16"
#endif
#if !defined(LEO_HAS_FF8)
#error "The AVX-512/GFNI T=16 kernel requires GF8"
#endif

namespace leopard { namespace backend {
namespace {

static const unsigned kSide = 16;
static const unsigned kShardBytes = 64;
static const uint8_t kZeroSkewLog = 255;

struct MultiplierMatrices
{
    uint64_t log7;
    uint64_t log17;
    uint64_t log28;
    uint64_t log34;
    uint64_t log51;
    uint64_t log85;
    uint64_t log102;
    uint64_t log111;
    uint64_t log131;
    uint64_t log153;
    uint64_t log183;
    uint64_t log187;
    uint64_t log219;
    uint64_t log222;
    uint64_t log224;
};

alignas(64) static MultiplierMatrices Matrices;

#if defined(__GNUC__) || defined(__clang__)
#define LEO2_T16_INLINE inline __attribute__((always_inline))
#if defined(__GNUC__) && !defined(__clang__) && defined(__ELF__)
#define LEO2_T16_KERNEL __attribute__((noinline, noipa, aligned(64), \
    section(".text.leo2_avx512_gfni_t16")))
#elif defined(__ELF__)
#define LEO2_T16_KERNEL __attribute__((noinline, aligned(64), \
    section(".text.leo2_avx512_gfni_t16")))
#else
#define LEO2_T16_KERNEL __attribute__((noinline, aligned(64)))
#endif
#else
#error "The AVX-512/GFNI T=16 object currently requires GCC or Clang"
#endif

static LEO2_T16_INLINE __m512i Broadcast(uint64_t matrix)
{
    return _mm512_set1_epi64(static_cast<long long>(matrix));
}

static LEO2_T16_INLINE __m512i Multiply(
    __m512i value, __m512i matrix)
{
    return _mm512_gf2p8affine_epi64_epi8(value, matrix, 0);
}

static LEO2_T16_INLINE void Inverse4(
    __m512i& x0,
    __m512i& x1,
    __m512i& x2,
    __m512i& x3,
    uint64_t matrix01_value,
    uint64_t matrix23_value,
    uint64_t matrix02_value)
{
    const __m512i matrix01 = Broadcast(matrix01_value);
    const __m512i matrix23 = Broadcast(matrix23_value);
    const __m512i matrix02 = Broadcast(matrix02_value);

    x1 = _mm512_xor_si512(x1, x0);
    x0 = _mm512_xor_si512(x0, Multiply(x1, matrix01));
    x3 = _mm512_xor_si512(x3, x2);
    x2 = _mm512_xor_si512(x2, Multiply(x3, matrix23));
    x2 = _mm512_xor_si512(x2, x0);
    x3 = _mm512_xor_si512(x3, x1);
    x0 = _mm512_xor_si512(x0, Multiply(x2, matrix02));
    x1 = _mm512_xor_si512(x1, Multiply(x3, matrix02));
}

static LEO2_T16_INLINE void Forward4(
    __m512i& x0,
    __m512i& x1,
    __m512i& x2,
    __m512i& x3,
    uint64_t matrix01_value,
    uint64_t matrix23_value,
    uint64_t matrix02_value)
{
    const __m512i matrix01 = Broadcast(matrix01_value);
    const __m512i matrix23 = Broadcast(matrix23_value);
    const __m512i matrix02 = Broadcast(matrix02_value);

    x0 = _mm512_xor_si512(x0, Multiply(x2, matrix02));
    x1 = _mm512_xor_si512(x1, Multiply(x3, matrix02));
    x2 = _mm512_xor_si512(x2, x0);
    x3 = _mm512_xor_si512(x3, x1);
    x0 = _mm512_xor_si512(x0, Multiply(x1, matrix01));
    x1 = _mm512_xor_si512(x1, x0);
    x2 = _mm512_xor_si512(x2, Multiply(x3, matrix23));
    x3 = _mm512_xor_si512(x3, x2);
}

static LEO2_T16_INLINE void Forward4MiddleOnly(
    __m512i& x0,
    __m512i& x1,
    __m512i& x2,
    __m512i& x3,
    uint64_t matrix23_value)
{
    const __m512i matrix23 = Broadcast(matrix23_value);

    x2 = _mm512_xor_si512(x2, x0);
    x3 = _mm512_xor_si512(x3, x1);
    x1 = _mm512_xor_si512(x1, x0);
    x2 = _mm512_xor_si512(x2, Multiply(x3, matrix23));
    x3 = _mm512_xor_si512(x3, x2);
}

/*
    Exact legacy-high K=R=T=16, B=64 transform.  One ZMM owns each row, so
    all sixteen rows remain register-resident from the first source load to
    the final recovery store.  The schedule is the generated AVX2 T16 circuit:
    four contiguous inverse radix-four groups, four fused outer inverse/forward
    groups, and four contiguous forward radix-four groups.
*/
static void LEO2_T16_KERNEL EncodeK16R16B64(
    const void* input_pointer,
    void* output_pointer)
{
    const uint8_t* const input = static_cast<const uint8_t*>(input_pointer);
    uint8_t* const output = static_cast<uint8_t*>(output_pointer);

    __m512i x0 = _mm512_loadu_si512(input + 0U * kShardBytes);
    __m512i x1 = _mm512_loadu_si512(input + 1U * kShardBytes);
    __m512i x2 = _mm512_loadu_si512(input + 2U * kShardBytes);
    __m512i x3 = _mm512_loadu_si512(input + 3U * kShardBytes);
    __m512i x4 = _mm512_loadu_si512(input + 4U * kShardBytes);
    __m512i x5 = _mm512_loadu_si512(input + 5U * kShardBytes);
    __m512i x6 = _mm512_loadu_si512(input + 6U * kShardBytes);
    __m512i x7 = _mm512_loadu_si512(input + 7U * kShardBytes);
    __m512i x8 = _mm512_loadu_si512(input + 8U * kShardBytes);
    __m512i x9 = _mm512_loadu_si512(input + 9U * kShardBytes);
    __m512i x10 = _mm512_loadu_si512(input + 10U * kShardBytes);
    __m512i x11 = _mm512_loadu_si512(input + 11U * kShardBytes);
    __m512i x12 = _mm512_loadu_si512(input + 12U * kShardBytes);
    __m512i x13 = _mm512_loadu_si512(input + 13U * kShardBytes);
    __m512i x14 = _mm512_loadu_si512(input + 14U * kShardBytes);
    __m512i x15 = _mm512_loadu_si512(input + 15U * kShardBytes);

#if defined(__GNUC__) || defined(__clang__)
    /*
        The input is immutable, so GCC otherwise rematerializes a few rows as
        memory-source XOR operands.  This empty definition barrier makes the
        sixteen loaded values the sole roots of the transform without emitting
        an instruction; the 32-register ZMM file still leaves ample temporaries.
    */
    __asm__ __volatile__(""
        : "+v"(x0), "+v"(x1), "+v"(x2), "+v"(x3),
          "+v"(x4), "+v"(x5), "+v"(x6), "+v"(x7));
    __asm__ __volatile__(""
        : "+v"(x8), "+v"(x9), "+v"(x10), "+v"(x11),
          "+v"(x12), "+v"(x13), "+v"(x14), "+v"(x15));
#endif

    Inverse4(x0, x1, x2, x3,
        Matrices.log219, Matrices.log7, Matrices.log153);
    Inverse4(x4, x5, x6, x7,
        Matrices.log111, Matrices.log28, Matrices.log102);
    Inverse4(x8, x9, x10, x11,
        Matrices.log183, Matrices.log224, Matrices.log51);
    Inverse4(x12, x13, x14, x15,
        Matrices.log131, Matrices.log222, Matrices.log187);

    Inverse4(x0, x4, x8, x12,
        Matrices.log17, Matrices.log34, Matrices.log85);
    Forward4MiddleOnly(x0, x4, x8, x12, Matrices.log85);
    Inverse4(x1, x5, x9, x13,
        Matrices.log17, Matrices.log34, Matrices.log85);
    Forward4MiddleOnly(x1, x5, x9, x13, Matrices.log85);
    Inverse4(x2, x6, x10, x14,
        Matrices.log17, Matrices.log34, Matrices.log85);
    Forward4MiddleOnly(x2, x6, x10, x14, Matrices.log85);
    Inverse4(x3, x7, x11, x15,
        Matrices.log17, Matrices.log34, Matrices.log85);
    Forward4MiddleOnly(x3, x7, x11, x15, Matrices.log85);

    Forward4MiddleOnly(x0, x1, x2, x3, Matrices.log85);
    Forward4(x4, x5, x6, x7,
        Matrices.log17, Matrices.log34, Matrices.log85);
    Forward4(x8, x9, x10, x11,
        Matrices.log153, Matrices.log102, Matrices.log17);
    Forward4(x12, x13, x14, x15,
        Matrices.log51, Matrices.log187, Matrices.log34);

    _mm512_storeu_si512(output + 0U * kShardBytes, x0);
    _mm512_storeu_si512(output + 1U * kShardBytes, x1);
    _mm512_storeu_si512(output + 2U * kShardBytes, x2);
    _mm512_storeu_si512(output + 3U * kShardBytes, x3);
    _mm512_storeu_si512(output + 4U * kShardBytes, x4);
    _mm512_storeu_si512(output + 5U * kShardBytes, x5);
    _mm512_storeu_si512(output + 6U * kShardBytes, x6);
    _mm512_storeu_si512(output + 7U * kShardBytes, x7);
    _mm512_storeu_si512(output + 8U * kShardBytes, x8);
    _mm512_storeu_si512(output + 9U * kShardBytes, x9);
    _mm512_storeu_si512(output + 10U * kShardBytes, x10);
    _mm512_storeu_si512(output + 11U * kShardBytes, x11);
    _mm512_storeu_si512(output + 12U * kShardBytes, x12);
    _mm512_storeu_si512(output + 13U * kShardBytes, x13);
    _mm512_storeu_si512(output + 14U * kShardBytes, x14);
    _mm512_storeu_si512(output + 15U * kShardBytes, x15);
}

static uint8_t ApplyMatrix(uint8_t value, uint64_t matrix)
{
    static const uint8_t kNibbleParity[16] = {
        0, 1, 1, 0, 1, 0, 0, 1,
        1, 0, 0, 1, 0, 1, 1, 0
    };
    uint8_t result = 0;
    for (unsigned output_bit = 0; output_bit < 8; ++output_bit)
    {
        const uint8_t row = static_cast<uint8_t>(
            matrix >> (8U * (7U - output_bit)));
        const uint8_t masked = static_cast<uint8_t>(row & value);
        const unsigned parity = static_cast<unsigned>(
            kNibbleParity[masked & 15U] ^ kNibbleParity[masked >> 4]);
        result |= static_cast<uint8_t>(parity << output_bit);
    }
    return result;
}

static bool BuildMatrix(
    FF8MultiplyLog multiply_log,
    uint8_t log,
    uint64_t& matrix)
{
    matrix = 0;
    for (unsigned output_bit = 0; output_bit < 8; ++output_bit)
    {
        uint8_t row = 0;
        for (unsigned input_bit = 0; input_bit < 8; ++input_bit)
        {
            const uint8_t product = multiply_log(
                static_cast<uint8_t>(1U << input_bit), log);
            if ((product & (1U << output_bit)) != 0)
                row |= static_cast<uint8_t>(1U << input_bit);
        }
        matrix |= static_cast<uint64_t>(row) <<
            (8U * (7U - output_bit));
    }
    for (unsigned value = 0; value < 256; ++value)
    {
        if (ApplyMatrix(static_cast<uint8_t>(value), matrix) !=
            multiply_log(static_cast<uint8_t>(value), log))
            return false;
    }
    return true;
}

static bool BuildMatrices(FF8MultiplyLog multiply_log)
{
    return multiply_log &&
        BuildMatrix(multiply_log, 7, Matrices.log7) &&
        BuildMatrix(multiply_log, 17, Matrices.log17) &&
        BuildMatrix(multiply_log, 28, Matrices.log28) &&
        BuildMatrix(multiply_log, 34, Matrices.log34) &&
        BuildMatrix(multiply_log, 51, Matrices.log51) &&
        BuildMatrix(multiply_log, 85, Matrices.log85) &&
        BuildMatrix(multiply_log, 102, Matrices.log102) &&
        BuildMatrix(multiply_log, 111, Matrices.log111) &&
        BuildMatrix(multiply_log, 131, Matrices.log131) &&
        BuildMatrix(multiply_log, 153, Matrices.log153) &&
        BuildMatrix(multiply_log, 183, Matrices.log183) &&
        BuildMatrix(multiply_log, 187, Matrices.log187) &&
        BuildMatrix(multiply_log, 219, Matrices.log219) &&
        BuildMatrix(multiply_log, 222, Matrices.log222) &&
        BuildMatrix(multiply_log, 224, Matrices.log224);
}

static void ScalarTransform(
    FF8MultiplyLog multiply_log,
    const uint8_t* skew_log_storage,
    const uint8_t* input,
    uint8_t* output,
    uint8_t state[kSide][kShardBytes])
{
    std::memcpy(state, input, sizeof(uint8_t) * kSide * kShardBytes);
    for (unsigned level = 0; level < 4; ++level)
    {
        const unsigned distance = 1U << level;
        const unsigned width = distance * 2U;
        for (unsigned row = 0; row < kSide; row += width)
        {
            const uint8_t log = skew_log_storage[kSide + row + distance];
            for (unsigned lane = 0; lane < distance; ++lane)
            {
                uint8_t* const x = state[row + lane];
                uint8_t* const y = state[row + lane + distance];
                for (unsigned column = 0; column < kShardBytes; ++column)
                {
                    y[column] ^= x[column];
                    if (log != kZeroSkewLog)
                        x[column] ^= multiply_log(y[column], log);
                }
            }
        }
    }
    for (unsigned level = 4; level-- > 0; )
    {
        const unsigned distance = 1U << level;
        const unsigned width = distance * 2U;
        for (unsigned row = 0; row < kSide; row += width)
        {
            const uint8_t log = skew_log_storage[row + distance];
            for (unsigned lane = 0; lane < distance; ++lane)
            {
                uint8_t* const x = state[row + lane];
                uint8_t* const y = state[row + lane + distance];
                for (unsigned column = 0; column < kShardBytes; ++column)
                {
                    if (log != kZeroSkewLog)
                        x[column] ^= multiply_log(y[column], log);
                    y[column] ^= x[column];
                }
            }
        }
    }
    std::memcpy(output, state, sizeof(uint8_t) * kSide * kShardBytes);
}

static bool KnownAnswerTest(
    FF8MultiplyLog multiply_log,
    const uint8_t* skew_log_storage)
{
    alignas(64) static uint8_t input[kSide * kShardBytes];
    alignas(64) static uint8_t expected[kSide * kShardBytes];
    alignas(64) static uint8_t actual[kSide * kShardBytes];
    alignas(64) static uint8_t scalar_state[kSide][kShardBytes];

    for (unsigned pattern = 0; pattern < 2; ++pattern)
    {
        uint64_t random = UINT64_C(0x74313667666e6935) + pattern;
        for (size_t i = 0; i < sizeof(input); ++i)
        {
            random += UINT64_C(0x9e3779b97f4a7c15);
            uint64_t value = random;
            value = (value ^ (value >> 30)) *
                UINT64_C(0xbf58476d1ce4e5b9);
            value = (value ^ (value >> 27)) *
                UINT64_C(0x94d049bb133111eb);
            value ^= value >> 31;
            input[i] = pattern == 0
                ? static_cast<uint8_t>(value)
                : static_cast<uint8_t>(i * 29U + (i >> 6) * 17U + 1U);
        }
        ScalarTransform(
            multiply_log, skew_log_storage, input, expected, scalar_state);
        EncodeK16R16B64(input, actual);
        if (std::memcmp(expected, actual, sizeof(expected)) != 0)
            return false;
    }
    return true;
}

} // namespace

FF8K16R16B64PackedKernel InitializeAVX512GFNIT16(
    FF8MultiplyLog multiply_log,
    const uint8_t* skew_log_storage)
{
    if (!multiply_log || !skew_log_storage ||
        !BuildMatrices(multiply_log) ||
        !KnownAnswerTest(multiply_log, skew_log_storage))
        return NULL;
    return &EncodeK16R16B64;
}

#undef LEO2_T16_KERNEL
#undef LEO2_T16_INLINE

}} // namespace leopard::backend
