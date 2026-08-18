/*
    Copyright (c) 2017 Christopher A. Taylor.  All rights reserved.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are met:

    * Redistributions of source code must retain the above copyright notice,
      this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright notice,
      this list of conditions and the following disclaimer in the documentation
      and/or other materials provided with the distribution.
    * Neither the name of Leopard-RS nor the names of its contributors may be
      used to endorse or promote products derived from this software without
      specific prior written permission.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
    AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
    IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
    ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
    LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
    CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
    SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
    INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
    CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
    ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
    POSSIBILITY OF SUCH DAMAGE.
*/

#include "Leopard2Backend.h"
#include "LeopardCommon.h"

#include <cstddef>
#include <cstdint>
#include <cstring>

#include <immintrin.h>

#if !defined(LEO2_HAVE_AVX512_GFNI_T128)
#error "This translation unit requires LEO2_HAVE_AVX512_GFNI_T128"
#endif
#if !defined(LEO_HAS_FF8)
#error "The AVX-512/GFNI T=128 kernel requires GF8"
#endif

namespace leopard { namespace backend {
namespace {

static const unsigned kOriginalCount = 65;
static const unsigned kRecoveryCount = 65;
static const unsigned kSide = 128;
static const unsigned kShardBytes = 64;
static const uint8_t kZeroSkewLog = 255;

alignas(64) static uint64_t MultiplierMatrixByLog[256];
alignas(64) static uint8_t SkewLogStorage[256];

#if defined(__GNUC__) || defined(__clang__)
#define LEO2_T128_INLINE inline __attribute__((always_inline))
#if defined(__GNUC__) && !defined(__clang__) && defined(__ELF__)
#define LEO2_T128_KERNEL __attribute__((noinline, noipa, aligned(64), \
    section(".text.leo2_avx512_gfni_t128")))
#elif defined(__ELF__)
#define LEO2_T128_KERNEL __attribute__((noinline, aligned(64), \
    section(".text.leo2_avx512_gfni_t128")))
#else
#define LEO2_T128_KERNEL __attribute__((noinline, aligned(64)))
#endif
#else
#error "The AVX-512/GFNI T=128 object currently requires GCC or Clang"
#endif

static LEO2_T128_INLINE __m512i LoadRow(
    const uint8_t* state, unsigned row)
{
    return _mm512_load_si512(reinterpret_cast<const __m512i*>(
        state + static_cast<size_t>(row) * kShardBytes));
}

static LEO2_T128_INLINE void StoreRow(
    uint8_t* state, unsigned row, __m512i value)
{
    _mm512_store_si512(reinterpret_cast<__m512i*>(
        state + static_cast<size_t>(row) * kShardBytes), value);
}

static LEO2_T128_INLINE __m512i Multiply(
    __m512i value, uint64_t matrix)
{
    return _mm512_gf2p8affine_epi64_epi8(
        value, _mm512_set1_epi64(static_cast<long long>(matrix)), 0);
}

static LEO2_T128_INLINE void Inverse4(
    uint8_t* state, unsigned row, unsigned distance,
    uint64_t matrix01, uint64_t matrix23, uint64_t matrix02)
{
    __m512i x0 = LoadRow(state, row);
    __m512i x1 = LoadRow(state, row + distance);
    __m512i x2 = LoadRow(state, row + distance * 2U);
    __m512i x3 = LoadRow(state, row + distance * 3U);

    x1 = _mm512_xor_si512(x1, x0);
    if (matrix01 != 0)
        x0 = _mm512_xor_si512(x0, Multiply(x1, matrix01));
    x3 = _mm512_xor_si512(x3, x2);
    if (matrix23 != 0)
        x2 = _mm512_xor_si512(x2, Multiply(x3, matrix23));
    x2 = _mm512_xor_si512(x2, x0);
    x3 = _mm512_xor_si512(x3, x1);
    if (matrix02 != 0)
    {
        x0 = _mm512_xor_si512(x0, Multiply(x2, matrix02));
        x1 = _mm512_xor_si512(x1, Multiply(x3, matrix02));
    }

    StoreRow(state, row, x0);
    StoreRow(state, row + distance, x1);
    StoreRow(state, row + distance * 2U, x2);
    StoreRow(state, row + distance * 3U, x3);
}

static LEO2_T128_INLINE void Forward4(
    uint8_t* state, unsigned row, unsigned distance,
    uint64_t matrix01, uint64_t matrix23, uint64_t matrix02)
{
    __m512i x0 = LoadRow(state, row);
    __m512i x1 = LoadRow(state, row + distance);
    __m512i x2 = LoadRow(state, row + distance * 2U);
    __m512i x3 = LoadRow(state, row + distance * 3U);

    if (matrix02 != 0)
    {
        x0 = _mm512_xor_si512(x0, Multiply(x2, matrix02));
        x1 = _mm512_xor_si512(x1, Multiply(x3, matrix02));
    }
    x2 = _mm512_xor_si512(x2, x0);
    x3 = _mm512_xor_si512(x3, x1);
    if (matrix01 != 0)
        x0 = _mm512_xor_si512(x0, Multiply(x1, matrix01));
    x1 = _mm512_xor_si512(x1, x0);
    if (matrix23 != 0)
        x2 = _mm512_xor_si512(x2, Multiply(x3, matrix23));
    x3 = _mm512_xor_si512(x3, x2);

    StoreRow(state, row, x0);
    StoreRow(state, row + distance, x1);
    StoreRow(state, row + distance * 2U, x2);
    StoreRow(state, row + distance * 3U, x3);
}

static LEO2_T128_INLINE void Inverse2(
    uint8_t* state, unsigned row, unsigned distance, uint64_t matrix)
{
    __m512i x = LoadRow(state, row);
    __m512i y = LoadRow(state, row + distance);
    y = _mm512_xor_si512(y, x);
    if (matrix != 0)
        x = _mm512_xor_si512(x, Multiply(y, matrix));
    StoreRow(state, row, x);
    StoreRow(state, row + distance, y);
}

static LEO2_T128_INLINE void Forward2(
    uint8_t* state, unsigned row, unsigned distance, uint64_t matrix)
{
    __m512i x = LoadRow(state, row);
    __m512i y = LoadRow(state, row + distance);
    if (matrix != 0)
        x = _mm512_xor_si512(x, Multiply(y, matrix));
    y = _mm512_xor_si512(y, x);
    StoreRow(state, row, x);
    StoreRow(state, row + distance, y);
}

static void LEO2_T128_KERNEL EncodeK65R65B64(
    const void* input_pointer,
    void* output_pointer,
    void* state_pointer)
{
    const uint8_t* const input = static_cast<const uint8_t*>(input_pointer);
    uint8_t* const output = static_cast<uint8_t*>(output_pointer);
    uint8_t* const state = static_cast<uint8_t*>(state_pointer);

    for (unsigned row = 0; row < kOriginalCount; ++row)
    {
        StoreRow(state, row, _mm512_loadu_si512(
            input + static_cast<size_t>(row) * kShardBytes));
    }
    const __m512i zero = _mm512_setzero_si512();
    for (unsigned row = kOriginalCount; row < kSide; ++row)
        StoreRow(state, row, zero);

    for (unsigned distance = 1; distance <= 16; distance *= 4)
    {
        const unsigned width = distance * 4U;
        for (unsigned row = 0; row < kOriginalCount; row += width)
        {
            const uint8_t* const skew = SkewLogStorage + 128U + row;
            const uint64_t matrix01 =
                MultiplierMatrixByLog[skew[distance]];
            const uint64_t matrix23 =
                MultiplierMatrixByLog[skew[distance * 3U]];
            const uint64_t matrix02 =
                MultiplierMatrixByLog[skew[distance * 2U]];
            for (unsigned lane = 0; lane < distance; ++lane)
            {
                Inverse4(state, row + lane, distance,
                    matrix01, matrix23, matrix02);
            }
        }
    }
    const uint64_t inverse_final =
        MultiplierMatrixByLog[SkewLogStorage[128U + 64U]];
    for (unsigned row = 0; row < 64; ++row)
        Inverse2(state, row, 64, inverse_final);

    for (unsigned distance = 32; distance >= 2; distance /= 4)
    {
        const unsigned width = distance * 4U;
        for (unsigned row = 0; row < kRecoveryCount; row += width)
        {
            const uint8_t* const skew = SkewLogStorage + row;
            const uint64_t matrix01 =
                MultiplierMatrixByLog[skew[distance]];
            const uint64_t matrix23 =
                MultiplierMatrixByLog[skew[distance * 3U]];
            const uint64_t matrix02 =
                MultiplierMatrixByLog[skew[distance * 2U]];
            for (unsigned lane = 0; lane < distance; ++lane)
            {
                Forward4(state, row + lane, distance,
                    matrix01, matrix23, matrix02);
            }
        }
    }
    for (unsigned row = 0; row < kRecoveryCount; row += 2)
    {
        Forward2(state, row, 1,
            MultiplierMatrixByLog[SkewLogStorage[row + 1U]]);
    }

    for (unsigned row = 0; row < kRecoveryCount; ++row)
    {
        _mm512_storeu_si512(
            output + static_cast<size_t>(row) * kShardBytes,
            LoadRow(state, row));
    }
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

static bool BuildMatrices(FF8MultiplyLog multiply_log)
{
    if (!multiply_log)
        return false;
    for (unsigned log = 0; log < 256; ++log)
    {
        if (log == kZeroSkewLog)
        {
            MultiplierMatrixByLog[log] = 0;
            continue;
        }
        uint64_t matrix = 0;
        for (unsigned output_bit = 0; output_bit < 8; ++output_bit)
        {
            uint8_t row = 0;
            for (unsigned input_bit = 0; input_bit < 8; ++input_bit)
            {
                const uint8_t product = multiply_log(
                    static_cast<uint8_t>(1U << input_bit),
                    static_cast<uint8_t>(log));
                if ((product & (1U << output_bit)) != 0)
                    row |= static_cast<uint8_t>(1U << input_bit);
            }
            matrix |= static_cast<uint64_t>(row) <<
                (8U * (7U - output_bit));
        }
        MultiplierMatrixByLog[log] = matrix;
        for (unsigned value = 0; value < 256; ++value)
        {
            if (ApplyMatrix(static_cast<uint8_t>(value), matrix) !=
                multiply_log(static_cast<uint8_t>(value),
                    static_cast<uint8_t>(log)))
                return false;
        }
    }
    return true;
}

static void ScalarTransform(
    FF8MultiplyLog multiply_log,
    const uint8_t* input,
    uint8_t* output,
    uint8_t state[128][64])
{
    std::memcpy(state, input, kOriginalCount * kShardBytes);
    std::memset(state[kOriginalCount], 0,
        (kSide - kOriginalCount) * kShardBytes);
    for (unsigned level = 0; level < 7; ++level)
    {
        const unsigned distance = 1U << level;
        const unsigned width = distance * 2U;
        for (unsigned row = 0; row < kSide; row += width)
        {
            const uint8_t log =
                SkewLogStorage[128U + row + distance];
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
    for (unsigned level = 7; level-- > 0; )
    {
        const unsigned distance = 1U << level;
        const unsigned width = distance * 2U;
        for (unsigned row = 0; row < kSide; row += width)
        {
            const uint8_t log = SkewLogStorage[row + distance];
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
    std::memcpy(output, state, kRecoveryCount * kShardBytes);
}

static bool KnownAnswerTest(FF8MultiplyLog multiply_log)
{
    alignas(64) static uint8_t input[kOriginalCount * kShardBytes];
    alignas(64) static uint8_t expected[kRecoveryCount * kShardBytes];
    alignas(64) static uint8_t actual[kRecoveryCount * kShardBytes];
    alignas(64) static uint8_t scalar_state[kSide][kShardBytes];
    alignas(64) static uint8_t vector_state[kSide * kShardBytes];

    for (unsigned pattern = 0; pattern < 2; ++pattern)
    {
        uint64_t random = UINT64_C(0x6b36357236356236) + pattern;
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
        ScalarTransform(multiply_log, input, expected, scalar_state);
        EncodeK65R65B64(input, actual, vector_state);
        if (std::memcmp(expected, actual, sizeof(expected)) != 0)
            return false;
    }
    return true;
}

} // namespace

FF8K65R65B64PackedKernel InitializeAVX512GFNIT128(
    FF8MultiplyLog multiply_log,
    const uint8_t* skew_log_storage)
{
    if (!multiply_log || !skew_log_storage)
        return NULL;
    std::memcpy(SkewLogStorage, skew_log_storage,
        sizeof(SkewLogStorage));
    if (!BuildMatrices(multiply_log) || !KnownAnswerTest(multiply_log))
        return NULL;
    return &EncodeK65R65B64;
}

#undef LEO2_T128_KERNEL
#undef LEO2_T128_INLINE

}} // namespace leopard::backend
