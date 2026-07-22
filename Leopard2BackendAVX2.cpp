/*
    Copyright (c) 2017 Christopher A. Taylor.  All rights reserved.
    See Leopard2Backend.h for the full BSD license notice.
*/

#include "Leopard2Backend.h"
#include "LeopardCommon.h"

#include <algorithm>
#include <cstring>
#include <immintrin.h>
#include <memory>
#include <new>

namespace leopard { namespace backend {

#if defined(LEO2_AVX512_VARIANT)
# define LEO2_AVX_BACKEND_KIND LEO2_BACKEND_AVX512
# define LEO2_AVX_BACKEND_NAME "avx512vl"
# define LEO2_AVX_INITIALIZER InitializeAVX512
# define LEO2_AVX_TABLE_STATE TestGetAVX512TableState
#else
# define LEO2_AVX_BACKEND_KIND LEO2_BACKEND_AVX2
# define LEO2_AVX_BACKEND_NAME "avx2"
# define LEO2_AVX_INITIALIZER InitializeAVX2
# define LEO2_AVX_TABLE_STATE TestGetAVX2TableState
#endif

#ifdef LEO_HAS_FF8
struct FF8NibbleTable
{
    uint8_t low[16];
    uint8_t high[16];
};
static const FF8NibbleTable* FF8Tables = NULL;
#endif

#ifdef LEO_HAS_FF16
struct FF16NibbleTable
{
    uint8_t low[4][16];
    uint8_t high[4][16];
};
static const FF16NibbleTable* FF16Tables = NULL;
#endif

#ifdef LEO_HAS_FF8
static uint8_t FF8Product(uint16_t log, uint8_t value)
{
    const FF8NibbleTable& table = FF8Tables[log];
    return static_cast<uint8_t>(
        table.low[value & 15U] ^ table.high[value >> 4]);
}
#endif

static __m256i BroadcastTable(const uint8_t table[16])
{
    return _mm256_broadcastsi128_si256(_mm_loadu_si128(
        reinterpret_cast<const __m128i*>(table)));
}

#ifdef LEO_HAS_FF8
template<bool Add>
static void AVX2FF8Operation(
    void* destination,
    const void* source,
    uint16_t multiplier_log,
    uint64_t byte_count)
{
    uint8_t* output = static_cast<uint8_t*>(destination);
    const uint8_t* input = static_cast<const uint8_t*>(source);
    const FF8NibbleTable& table = FF8Tables[multiplier_log];
    const __m256i low_table = BroadcastTable(table.low);
    const __m256i high_table = BroadcastTable(table.high);
    const __m256i nibble_mask = _mm256_set1_epi8(15);
    while (byte_count >= 32)
    {
        const __m256i data = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(input));
        const __m256i low = _mm256_shuffle_epi8(
            low_table, _mm256_and_si256(data, nibble_mask));
        const __m256i high = _mm256_shuffle_epi8(high_table,
            _mm256_and_si256(_mm256_srli_epi64(data, 4), nibble_mask));
        __m256i product = _mm256_xor_si256(low, high);
        if (Add)
            product = _mm256_xor_si256(product, _mm256_loadu_si256(
                reinterpret_cast<const __m256i*>(output)));
        _mm256_storeu_si256(reinterpret_cast<__m256i*>(output), product);
        input += 32;
        output += 32;
        byte_count -= 32;
    }
    while (byte_count-- != 0)
    {
        const uint8_t product = FF8Product(multiplier_log, *input++);
        if (Add)
            *output++ ^= product;
        else
            *output++ = product;
    }
}

static void AVX2FF8Multiply(
    void* destination, const void* source,
    uint16_t multiplier_log, uint64_t byte_count)
{
    AVX2FF8Operation<false>(destination, source, multiplier_log, byte_count);
}

static void AVX2FF8MultiplyAdd(
    void* destination, const void* source,
    uint16_t multiplier_log, uint64_t byte_count)
{
    AVX2FF8Operation<true>(destination, source, multiplier_log, byte_count);
}

static __m256i AVX2FF8ProductVector(
    __m256i data,
    __m256i low_table,
    __m256i high_table)
{
    const __m256i nibble_mask = _mm256_set1_epi8(15);
    const __m256i low = _mm256_shuffle_epi8(
        low_table, _mm256_and_si256(data, nibble_mask));
    const __m256i high = _mm256_shuffle_epi8(high_table,
        _mm256_and_si256(_mm256_srli_epi64(data, 4), nibble_mask));
    return _mm256_xor_si256(low, high);
}

static __m256i AVX2FF8ApplyWeight(__m256i data, uint16_t weight_log)
{
    static const uint16_t kModulus = 255;
    if (weight_log == 0 || weight_log == kModulus)
        return data;
    const FF8NibbleTable& table = FF8Tables[weight_log];
    return AVX2FF8ProductVector(
        data, BroadcastTable(table.low), BroadcastTable(table.high));
}

template<bool Inverse>
static void AVX2FF8Butterfly2(
    void* x_pointer,
    void* y_pointer,
    uint16_t multiplier_log,
    uint64_t byte_count)
{
    uint8_t* x = static_cast<uint8_t*>(x_pointer);
    uint8_t* y = static_cast<uint8_t*>(y_pointer);
    const FF8NibbleTable& table = FF8Tables[multiplier_log];
    const __m256i low_table = BroadcastTable(table.low);
    const __m256i high_table = BroadcastTable(table.high);
    while (byte_count >= 32)
    {
        __m256i x_value = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(x));
        __m256i y_value = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(y));
        if (Inverse)
        {
            y_value = _mm256_xor_si256(y_value, x_value);
            x_value = _mm256_xor_si256(x_value,
                AVX2FF8ProductVector(y_value, low_table, high_table));
        }
        else
        {
            x_value = _mm256_xor_si256(x_value,
                AVX2FF8ProductVector(y_value, low_table, high_table));
            y_value = _mm256_xor_si256(y_value, x_value);
        }
        _mm256_storeu_si256(reinterpret_cast<__m256i*>(x), x_value);
        _mm256_storeu_si256(reinterpret_cast<__m256i*>(y), y_value);
        x += 32;
        y += 32;
        byte_count -= 32;
    }
    while (byte_count-- != 0)
    {
        if (Inverse)
        {
            *y ^= *x;
            *x ^= FF8Product(multiplier_log, *y);
        }
        else
        {
            *x ^= FF8Product(multiplier_log, *y);
            *y ^= *x;
        }
        ++x;
        ++y;
    }
}

static void AVX2FF8IFFTButterfly2(
    void* x, void* y, uint16_t multiplier_log, uint64_t byte_count)
{
    AVX2FF8Butterfly2<true>(x, y, multiplier_log, byte_count);
}

static void AVX2FF8FFTButterfly2(
    void* x, void* y, uint16_t multiplier_log, uint64_t byte_count)
{
    AVX2FF8Butterfly2<false>(x, y, multiplier_log, byte_count);
}

static void AVX2FF8FFTButterfly2Out(
    const void* x_input_pointer,
    const void* y_input_pointer,
    void* x_output_pointer,
    void* y_output_pointer,
    uint16_t multiplier_log,
    uint64_t byte_count)
{
    static const uint16_t kZeroSkew = 255;
    const uint8_t* x_input = static_cast<const uint8_t*>(x_input_pointer);
    const uint8_t* y_input = static_cast<const uint8_t*>(y_input_pointer);
    uint8_t* x_output = static_cast<uint8_t*>(x_output_pointer);
    uint8_t* y_output = static_cast<uint8_t*>(y_output_pointer);
    __m256i low_table = _mm256_setzero_si256();
    __m256i high_table = _mm256_setzero_si256();
    if (multiplier_log != kZeroSkew)
    {
        low_table = BroadcastTable(FF8Tables[multiplier_log].low);
        high_table = BroadcastTable(FF8Tables[multiplier_log].high);
    }
    while (byte_count >= 32)
    {
        __m256i x_value = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(x_input));
        __m256i y_value = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(y_input));
        if (multiplier_log != kZeroSkew)
            x_value = _mm256_xor_si256(x_value,
                AVX2FF8ProductVector(y_value, low_table, high_table));
        y_value = _mm256_xor_si256(y_value, x_value);
        _mm256_storeu_si256(reinterpret_cast<__m256i*>(x_output), x_value);
        _mm256_storeu_si256(reinterpret_cast<__m256i*>(y_output), y_value);
        x_input += 32;
        y_input += 32;
        x_output += 32;
        y_output += 32;
        byte_count -= 32;
    }
    while (byte_count-- != 0)
    {
        uint8_t x_value = *x_input++;
        uint8_t y_value = *y_input++;
        if (multiplier_log != kZeroSkew)
            x_value ^= FF8Product(multiplier_log, y_value);
        y_value ^= x_value;
        *x_output++ = x_value;
        *y_output++ = y_value;
    }
}

#if !defined(LEO2_AVX512_VARIANT)
static void AVX2FF8IFFTButterfly2Out(
    const void* x_input_pointer,
    const void* y_input_pointer,
    void* x_output_pointer,
    void* y_output_pointer,
    uint16_t multiplier_log,
    uint64_t byte_count)
{
    static const uint16_t kZeroSkew = 255;
    const uint8_t* x_input = static_cast<const uint8_t*>(x_input_pointer);
    const uint8_t* y_input = static_cast<const uint8_t*>(y_input_pointer);
    uint8_t* x_output = static_cast<uint8_t*>(x_output_pointer);
    uint8_t* y_output = static_cast<uint8_t*>(y_output_pointer);
    __m256i low_table = _mm256_setzero_si256();
    __m256i high_table = _mm256_setzero_si256();
    if (multiplier_log != kZeroSkew)
    {
        low_table = BroadcastTable(FF8Tables[multiplier_log].low);
        high_table = BroadcastTable(FF8Tables[multiplier_log].high);
    }
    while (byte_count >= 32)
    {
        const __m256i x_original = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(x_input));
        const __m256i y_value = _mm256_xor_si256(x_original,
            _mm256_loadu_si256(reinterpret_cast<const __m256i*>(y_input)));
        __m256i x_value = x_original;
        if (multiplier_log != kZeroSkew)
            x_value = _mm256_xor_si256(x_value,
                AVX2FF8ProductVector(y_value, low_table, high_table));
        _mm256_storeu_si256(reinterpret_cast<__m256i*>(x_output), x_value);
        _mm256_storeu_si256(reinterpret_cast<__m256i*>(y_output), y_value);
        x_input += 32;
        y_input += 32;
        x_output += 32;
        y_output += 32;
        byte_count -= 32;
    }
    while (byte_count-- != 0)
    {
        const uint8_t y_value = static_cast<uint8_t>(*y_input ^ *x_input);
        uint8_t x_value = *x_input;
        if (multiplier_log != kZeroSkew)
            x_value ^= FF8Product(multiplier_log, y_value);
        *x_output++ = x_value;
        *y_output++ = y_value;
        ++x_input;
        ++y_input;
    }
}
#endif

static void AVX2FF8IFFTButterfly2Xor(
    const void* x_input_pointer,
    const void* y_input_pointer,
    void* x_output_pointer,
    void* y_output_pointer,
    uint16_t multiplier_log,
    uint64_t byte_count)
{
    const uint8_t* x_input = static_cast<const uint8_t*>(x_input_pointer);
    const uint8_t* y_input = static_cast<const uint8_t*>(y_input_pointer);
    uint8_t* x_output = static_cast<uint8_t*>(x_output_pointer);
    uint8_t* y_output = static_cast<uint8_t*>(y_output_pointer);
    const FF8NibbleTable& table = FF8Tables[multiplier_log];
    const __m256i low_table = BroadcastTable(table.low);
    const __m256i high_table = BroadcastTable(table.high);
    while (byte_count >= 32)
    {
        const __m256i x_original = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(x_input));
        const __m256i y_value = _mm256_xor_si256(x_original,
            _mm256_loadu_si256(reinterpret_cast<const __m256i*>(y_input)));
        const __m256i x_value = _mm256_xor_si256(x_original,
            AVX2FF8ProductVector(y_value, low_table, high_table));
        const __m256i x_result = _mm256_xor_si256(x_value,
            _mm256_loadu_si256(reinterpret_cast<const __m256i*>(x_output)));
        const __m256i y_result = _mm256_xor_si256(y_value,
            _mm256_loadu_si256(reinterpret_cast<const __m256i*>(y_output)));
        _mm256_storeu_si256(reinterpret_cast<__m256i*>(x_output), x_result);
        _mm256_storeu_si256(reinterpret_cast<__m256i*>(y_output), y_result);
        x_input += 32;
        y_input += 32;
        x_output += 32;
        y_output += 32;
        byte_count -= 32;
    }
    while (byte_count-- != 0)
    {
        const uint8_t y_value = static_cast<uint8_t>(*y_input ^ *x_input);
        const uint8_t x_value = static_cast<uint8_t>(
            *x_input ^ FF8Product(multiplier_log, y_value));
        *x_output ^= x_value;
        *y_output ^= y_value;
        ++x_input;
        ++y_input;
        ++x_output;
        ++y_output;
    }
}

#endif // LEO_HAS_FF8

#ifdef LEO_HAS_FF16
static uint16_t FF16Product(uint16_t log, uint16_t value)
{
    const FF16NibbleTable& table = FF16Tables[log];
    const unsigned n0 = value & 15U;
    const unsigned n1 = (value >> 4) & 15U;
    const unsigned n2 = (value >> 8) & 15U;
    const unsigned n3 = value >> 12;
    const uint8_t low = static_cast<uint8_t>(
        table.low[0][n0] ^ table.low[1][n1] ^
        table.low[2][n2] ^ table.low[3][n3]);
    const uint8_t high = static_cast<uint8_t>(
        table.high[0][n0] ^ table.high[1][n1] ^
        table.high[2][n2] ^ table.high[3][n3]);
    return static_cast<uint16_t>(low | (static_cast<unsigned>(high) << 8));
}

static void AVX2FF16ProductVectors(
    __m256i low_data,
    __m256i high_data,
    const __m256i low_tables[4],
    const __m256i high_tables[4],
    __m256i& product_low,
    __m256i& product_high)
{
    const __m256i mask = _mm256_set1_epi8(15);
    const __m256i nibbles[4] = {
        _mm256_and_si256(low_data, mask),
        _mm256_and_si256(_mm256_srli_epi64(low_data, 4), mask),
        _mm256_and_si256(high_data, mask),
        _mm256_and_si256(_mm256_srli_epi64(high_data, 4), mask)
    };
    product_low = _mm256_setzero_si256();
    product_high = _mm256_setzero_si256();
    for (unsigned i = 0; i < 4; ++i)
    {
        product_low = _mm256_xor_si256(product_low,
            _mm256_shuffle_epi8(low_tables[i], nibbles[i]));
        product_high = _mm256_xor_si256(product_high,
            _mm256_shuffle_epi8(high_tables[i], nibbles[i]));
    }
}

template<bool Add>
static void AVX2FF16Operation(
    void* destination,
    const void* source,
    uint16_t multiplier_log,
    uint64_t byte_count)
{
    uint8_t* output = static_cast<uint8_t*>(destination);
    const uint8_t* input = static_cast<const uint8_t*>(source);
    const FF16NibbleTable& table = FF16Tables[multiplier_log];
    const __m256i low_tables[4] = {
        BroadcastTable(table.low[0]), BroadcastTable(table.low[1]),
        BroadcastTable(table.low[2]), BroadcastTable(table.low[3])
    };
    const __m256i high_tables[4] = {
        BroadcastTable(table.high[0]), BroadcastTable(table.high[1]),
        BroadcastTable(table.high[2]), BroadcastTable(table.high[3])
    };
    const __m256i mask = _mm256_set1_epi8(15);
    uint64_t offset = 0;
    while (byte_count - offset >= 64)
    {
        const __m256i low_data = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(input + offset));
        const __m256i high_data = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(input + offset + 32));
        const __m256i nibbles[4] = {
            _mm256_and_si256(low_data, mask),
            _mm256_and_si256(_mm256_srli_epi64(low_data, 4), mask),
            _mm256_and_si256(high_data, mask),
            _mm256_and_si256(_mm256_srli_epi64(high_data, 4), mask)
        };
        __m256i product_low = _mm256_setzero_si256();
        __m256i product_high = _mm256_setzero_si256();
        for (unsigned i = 0; i < 4; ++i)
        {
            product_low = _mm256_xor_si256(product_low,
                _mm256_shuffle_epi8(low_tables[i], nibbles[i]));
            product_high = _mm256_xor_si256(product_high,
                _mm256_shuffle_epi8(high_tables[i], nibbles[i]));
        }
        if (Add)
        {
            product_low = _mm256_xor_si256(product_low,
                _mm256_loadu_si256(reinterpret_cast<const __m256i*>(
                    output + offset)));
            product_high = _mm256_xor_si256(product_high,
                _mm256_loadu_si256(reinterpret_cast<const __m256i*>(
                    output + offset + 32)));
        }
        _mm256_storeu_si256(
            reinterpret_cast<__m256i*>(output + offset), product_low);
        _mm256_storeu_si256(
            reinterpret_cast<__m256i*>(output + offset + 32), product_high);
        offset += 64;
    }

    const uint64_t symbols = (byte_count - offset) / 2;
    for (uint64_t i = 0; i < symbols; ++i)
    {
        const uint16_t value = static_cast<uint16_t>(input[offset + i] |
            (static_cast<unsigned>(input[offset + symbols + i]) << 8));
        const uint16_t product = FF16Product(multiplier_log, value);
        if (Add)
        {
            output[offset + i] ^= static_cast<uint8_t>(product);
            output[offset + symbols + i] ^=
                static_cast<uint8_t>(product >> 8);
        }
        else
        {
            output[offset + i] = static_cast<uint8_t>(product);
            output[offset + symbols + i] = static_cast<uint8_t>(product >> 8);
        }
    }
}

static void AVX2FF16Multiply(
    void* destination, const void* source,
    uint16_t multiplier_log, uint64_t byte_count)
{
    AVX2FF16Operation<false>(destination, source, multiplier_log, byte_count);
}

static void AVX2FF16MultiplyAdd(
    void* destination, const void* source,
    uint16_t multiplier_log, uint64_t byte_count)
{
    AVX2FF16Operation<true>(destination, source, multiplier_log, byte_count);
}

template<bool Inverse>
static LEO_FORCE_INLINE void AVX2FF16Butterfly2Prepared(
    void* x_pointer,
    void* y_pointer,
    uint16_t multiplier_log,
    const __m256i low_tables[4],
    const __m256i high_tables[4],
    uint64_t byte_count)
{
    uint8_t* x = static_cast<uint8_t*>(x_pointer);
    uint8_t* y = static_cast<uint8_t*>(y_pointer);
    uint64_t offset = 0;
    while (byte_count - offset >= 64)
    {
        __m256i x_low = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(x + offset));
        __m256i x_high = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(x + offset + 32));
        __m256i y_low = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(y + offset));
        __m256i y_high = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(y + offset + 32));
        if (Inverse)
        {
            y_low = _mm256_xor_si256(y_low, x_low);
            y_high = _mm256_xor_si256(y_high, x_high);
        }
        __m256i product_low;
        __m256i product_high;
        AVX2FF16ProductVectors(y_low, y_high, low_tables, high_tables,
            product_low, product_high);
        x_low = _mm256_xor_si256(x_low, product_low);
        x_high = _mm256_xor_si256(x_high, product_high);
        if (!Inverse)
        {
            y_low = _mm256_xor_si256(y_low, x_low);
            y_high = _mm256_xor_si256(y_high, x_high);
        }
        _mm256_storeu_si256(reinterpret_cast<__m256i*>(x + offset), x_low);
        _mm256_storeu_si256(
            reinterpret_cast<__m256i*>(x + offset + 32), x_high);
        _mm256_storeu_si256(reinterpret_cast<__m256i*>(y + offset), y_low);
        _mm256_storeu_si256(
            reinterpret_cast<__m256i*>(y + offset + 32), y_high);
        offset += 64;
    }
    const uint64_t symbols = (byte_count - offset) / 2;
    for (uint64_t i = 0; i < symbols; ++i)
    {
        uint16_t x_value = static_cast<uint16_t>(x[offset + i] |
            (static_cast<unsigned>(x[offset + symbols + i]) << 8));
        uint16_t y_value = static_cast<uint16_t>(y[offset + i] |
            (static_cast<unsigned>(y[offset + symbols + i]) << 8));
        if (Inverse)
        {
            y_value ^= x_value;
            x_value ^= FF16Product(multiplier_log, y_value);
        }
        else
        {
            x_value ^= FF16Product(multiplier_log, y_value);
            y_value ^= x_value;
        }
        x[offset + i] = static_cast<uint8_t>(x_value);
        x[offset + symbols + i] = static_cast<uint8_t>(x_value >> 8);
        y[offset + i] = static_cast<uint8_t>(y_value);
        y[offset + symbols + i] = static_cast<uint8_t>(y_value >> 8);
    }
}

template<bool Inverse>
static void AVX2FF16Butterfly2(
    void* x_pointer,
    void* y_pointer,
    uint16_t multiplier_log,
    uint64_t byte_count)
{
    const FF16NibbleTable& table = FF16Tables[multiplier_log];
    const __m256i low_tables[4] = {
        BroadcastTable(table.low[0]), BroadcastTable(table.low[1]),
        BroadcastTable(table.low[2]), BroadcastTable(table.low[3])
    };
    const __m256i high_tables[4] = {
        BroadcastTable(table.high[0]), BroadcastTable(table.high[1]),
        BroadcastTable(table.high[2]), BroadcastTable(table.high[3])
    };
    AVX2FF16Butterfly2Prepared<Inverse>(
        x_pointer, y_pointer, multiplier_log,
        low_tables, high_tables, byte_count);
}

static void AVX2FF16IFFTButterfly2(
    void* x, void* y, uint16_t multiplier_log, uint64_t byte_count)
{
    AVX2FF16Butterfly2<true>(x, y, multiplier_log, byte_count);
}

static void AVX2FF16FFTButterfly2(
    void* x, void* y, uint16_t multiplier_log, uint64_t byte_count)
{
    AVX2FF16Butterfly2<false>(x, y, multiplier_log, byte_count);
}

static void AVX2FF16FFTButterfly2Out(
    const void* x_input_pointer,
    const void* y_input_pointer,
    void* x_output_pointer,
    void* y_output_pointer,
    uint16_t multiplier_log,
    uint64_t byte_count)
{
    static const uint16_t kZeroSkew = 65535;
    const uint8_t* x_input = static_cast<const uint8_t*>(x_input_pointer);
    const uint8_t* y_input = static_cast<const uint8_t*>(y_input_pointer);
    uint8_t* x_output = static_cast<uint8_t*>(x_output_pointer);
    uint8_t* y_output = static_cast<uint8_t*>(y_output_pointer);
    __m256i low_tables[4];
    __m256i high_tables[4];
    if (multiplier_log != kZeroSkew)
    {
        const FF16NibbleTable& table = FF16Tables[multiplier_log];
        for (unsigned i = 0; i < 4; ++i)
        {
            low_tables[i] = BroadcastTable(table.low[i]);
            high_tables[i] = BroadcastTable(table.high[i]);
        }
    }
    uint64_t offset = 0;
    while (byte_count - offset >= 64)
    {
        __m256i x_low = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(x_input + offset));
        __m256i x_high = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(x_input + offset + 32));
        __m256i y_low = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(y_input + offset));
        __m256i y_high = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(y_input + offset + 32));
        if (multiplier_log != kZeroSkew)
        {
            __m256i product_low;
            __m256i product_high;
            AVX2FF16ProductVectors(y_low, y_high, low_tables, high_tables,
                product_low, product_high);
            x_low = _mm256_xor_si256(x_low, product_low);
            x_high = _mm256_xor_si256(x_high, product_high);
        }
        y_low = _mm256_xor_si256(y_low, x_low);
        y_high = _mm256_xor_si256(y_high, x_high);
        _mm256_storeu_si256(
            reinterpret_cast<__m256i*>(x_output + offset), x_low);
        _mm256_storeu_si256(
            reinterpret_cast<__m256i*>(x_output + offset + 32), x_high);
        _mm256_storeu_si256(
            reinterpret_cast<__m256i*>(y_output + offset), y_low);
        _mm256_storeu_si256(
            reinterpret_cast<__m256i*>(y_output + offset + 32), y_high);
        offset += 64;
    }
    const uint64_t symbols = (byte_count - offset) / 2;
    for (uint64_t i = 0; i < symbols; ++i)
    {
        uint16_t x_value = static_cast<uint16_t>(x_input[offset + i] |
            (static_cast<unsigned>(x_input[offset + symbols + i]) << 8));
        uint16_t y_value = static_cast<uint16_t>(y_input[offset + i] |
            (static_cast<unsigned>(y_input[offset + symbols + i]) << 8));
        if (multiplier_log != kZeroSkew)
            x_value ^= FF16Product(multiplier_log, y_value);
        y_value ^= x_value;
        x_output[offset + i] = static_cast<uint8_t>(x_value);
        x_output[offset + symbols + i] = static_cast<uint8_t>(x_value >> 8);
        y_output[offset + i] = static_cast<uint8_t>(y_value);
        y_output[offset + symbols + i] = static_cast<uint8_t>(y_value >> 8);
    }
}

static void AVX2FF16IFFTButterfly2Xor(
    const void* x_input_pointer,
    const void* y_input_pointer,
    void* x_output_pointer,
    void* y_output_pointer,
    uint16_t multiplier_log,
    uint64_t byte_count)
{
    const uint8_t* x_input = static_cast<const uint8_t*>(x_input_pointer);
    const uint8_t* y_input = static_cast<const uint8_t*>(y_input_pointer);
    uint8_t* x_output = static_cast<uint8_t*>(x_output_pointer);
    uint8_t* y_output = static_cast<uint8_t*>(y_output_pointer);
    const FF16NibbleTable& table = FF16Tables[multiplier_log];
    const __m256i low_tables[4] = {
        BroadcastTable(table.low[0]), BroadcastTable(table.low[1]),
        BroadcastTable(table.low[2]), BroadcastTable(table.low[3])
    };
    const __m256i high_tables[4] = {
        BroadcastTable(table.high[0]), BroadcastTable(table.high[1]),
        BroadcastTable(table.high[2]), BroadcastTable(table.high[3])
    };
    uint64_t offset = 0;
    while (byte_count - offset >= 64)
    {
        __m256i x_low = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(x_input + offset));
        __m256i x_high = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(x_input + offset + 32));
        __m256i y_low = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(y_input + offset));
        __m256i y_high = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(y_input + offset + 32));
        y_low = _mm256_xor_si256(y_low, x_low);
        y_high = _mm256_xor_si256(y_high, x_high);
        __m256i product_low;
        __m256i product_high;
        AVX2FF16ProductVectors(y_low, y_high, low_tables, high_tables,
            product_low, product_high);
        x_low = _mm256_xor_si256(x_low, product_low);
        x_high = _mm256_xor_si256(x_high, product_high);
        x_low = _mm256_xor_si256(x_low, _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(x_output + offset)));
        x_high = _mm256_xor_si256(x_high, _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(x_output + offset + 32)));
        y_low = _mm256_xor_si256(y_low, _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(y_output + offset)));
        y_high = _mm256_xor_si256(y_high, _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(y_output + offset + 32)));
        _mm256_storeu_si256(
            reinterpret_cast<__m256i*>(x_output + offset), x_low);
        _mm256_storeu_si256(
            reinterpret_cast<__m256i*>(x_output + offset + 32), x_high);
        _mm256_storeu_si256(
            reinterpret_cast<__m256i*>(y_output + offset), y_low);
        _mm256_storeu_si256(
            reinterpret_cast<__m256i*>(y_output + offset + 32), y_high);
        offset += 64;
    }
    const uint64_t symbols = (byte_count - offset) / 2;
    for (uint64_t i = 0; i < symbols; ++i)
    {
        uint16_t x_value = static_cast<uint16_t>(x_input[offset + i] |
            (static_cast<unsigned>(x_input[offset + symbols + i]) << 8));
        uint16_t y_value = static_cast<uint16_t>(y_input[offset + i] |
            (static_cast<unsigned>(y_input[offset + symbols + i]) << 8));
        y_value ^= x_value;
        x_value ^= FF16Product(multiplier_log, y_value);
        x_output[offset + i] ^= static_cast<uint8_t>(x_value);
        x_output[offset + symbols + i] ^=
            static_cast<uint8_t>(x_value >> 8);
        y_output[offset + i] ^= static_cast<uint8_t>(y_value);
        y_output[offset + symbols + i] ^=
            static_cast<uint8_t>(y_value >> 8);
    }
}

#endif // LEO_HAS_FF16

static void AVX2XorMemory(
    void* destination, const void* source, uint64_t byte_count)
{
    uint8_t* output = static_cast<uint8_t*>(destination);
    const uint8_t* input = static_cast<const uint8_t*>(source);
    while (byte_count >= 32)
    {
        const __m256i result = _mm256_xor_si256(
            _mm256_loadu_si256(reinterpret_cast<const __m256i*>(output)),
            _mm256_loadu_si256(reinterpret_cast<const __m256i*>(input)));
        _mm256_storeu_si256(reinterpret_cast<__m256i*>(output), result);
        output += 32;
        input += 32;
        byte_count -= 32;
    }
    while (byte_count-- != 0)
        *output++ ^= *input++;
}

static void AVX2XorMemory2To1(
    void* destination,
    const void* source0,
    const void* source1,
    uint64_t byte_count)
{
    uint8_t* output = static_cast<uint8_t*>(destination);
    const uint8_t* input0 = static_cast<const uint8_t*>(source0);
    const uint8_t* input1 = static_cast<const uint8_t*>(source1);
    while (byte_count >= 128)
    {
        for (unsigned offset = 0; offset < 128; offset += 32)
        {
            __m256i result = _mm256_xor_si256(
                _mm256_loadu_si256(reinterpret_cast<const __m256i*>(
                    output + offset)),
                _mm256_loadu_si256(reinterpret_cast<const __m256i*>(
                    input0 + offset)));
            result = _mm256_xor_si256(result,
                _mm256_loadu_si256(reinterpret_cast<const __m256i*>(
                    input1 + offset)));
            _mm256_storeu_si256(
                reinterpret_cast<__m256i*>(output + offset), result);
        }
        output += 128;
        input0 += 128;
        input1 += 128;
        byte_count -= 128;
    }
    while (byte_count >= 32)
    {
        __m256i result = _mm256_xor_si256(
            _mm256_loadu_si256(reinterpret_cast<const __m256i*>(output)),
            _mm256_loadu_si256(reinterpret_cast<const __m256i*>(input0)));
        result = _mm256_xor_si256(result,
            _mm256_loadu_si256(reinterpret_cast<const __m256i*>(input1)));
        _mm256_storeu_si256(reinterpret_cast<__m256i*>(output), result);
        output += 32;
        input0 += 32;
        input1 += 32;
        byte_count -= 32;
    }
    while (byte_count-- != 0)
        *output++ ^= *input0++ ^ *input1++;
}

static void AVX2XorMemorySources(
    void* destination,
    const void* initial_source,
    const void* const* sources,
    uint32_t source_count,
    uint64_t byte_count)
{
    std::memcpy(destination, initial_source, static_cast<size_t>(byte_count));
    const void* waiting[6];
    unsigned waiting_count = 0;
    for (uint32_t i = 0; i < source_count; ++i)
    {
        if (!sources[i])
            continue;
        waiting[waiting_count++] = sources[i];
        if (waiting_count == 6)
        {
            uint8_t* output = static_cast<uint8_t*>(destination);
            const uint8_t* input[6] = {
                static_cast<const uint8_t*>(waiting[0]),
                static_cast<const uint8_t*>(waiting[1]),
                static_cast<const uint8_t*>(waiting[2]),
                static_cast<const uint8_t*>(waiting[3]),
                static_cast<const uint8_t*>(waiting[4]),
                static_cast<const uint8_t*>(waiting[5])
            };
            uint64_t offset = 0;
            while (byte_count - offset >= 32)
            {
                __m256i result = _mm256_loadu_si256(
                    reinterpret_cast<const __m256i*>(output + offset));
                for (unsigned lane = 0; lane < 6; ++lane)
                    result = _mm256_xor_si256(result, _mm256_loadu_si256(
                        reinterpret_cast<const __m256i*>(
                            input[lane] + offset)));
                _mm256_storeu_si256(
                    reinterpret_cast<__m256i*>(output + offset), result);
                offset += 32;
            }
            while (offset < byte_count)
            {
                output[offset] ^= input[0][offset] ^ input[1][offset] ^
                    input[2][offset] ^ input[3][offset] ^ input[4][offset] ^
                    input[5][offset];
                ++offset;
            }
            waiting_count = 0;
        }
    }
    while (waiting_count >= 2)
    {
        AVX2XorMemory2To1(
            destination, waiting[0], waiting[1], byte_count);
        waiting_count -= 2;
        for (unsigned i = 0; i < waiting_count; ++i)
            waiting[i] = waiting[i + 2];
    }
    if (waiting_count != 0)
        AVX2XorMemory(destination, waiting[0], byte_count);
}

static void AVX2XorMemory4(
    void* destination0, const void* source0,
    void* destination1, const void* source1,
    void* destination2, const void* source2,
    void* destination3, const void* source3,
    uint64_t byte_count)
{
    uint8_t* output0 = static_cast<uint8_t*>(destination0);
    uint8_t* output1 = static_cast<uint8_t*>(destination1);
    uint8_t* output2 = static_cast<uint8_t*>(destination2);
    uint8_t* output3 = static_cast<uint8_t*>(destination3);
    const uint8_t* input0 = static_cast<const uint8_t*>(source0);
    const uint8_t* input1 = static_cast<const uint8_t*>(source1);
    const uint8_t* input2 = static_cast<const uint8_t*>(source2);
    const uint8_t* input3 = static_cast<const uint8_t*>(source3);
    while (byte_count >= 32)
    {
        const __m256i result0 = _mm256_xor_si256(
            _mm256_loadu_si256(reinterpret_cast<const __m256i*>(output0)),
            _mm256_loadu_si256(reinterpret_cast<const __m256i*>(input0)));
        const __m256i result1 = _mm256_xor_si256(
            _mm256_loadu_si256(reinterpret_cast<const __m256i*>(output1)),
            _mm256_loadu_si256(reinterpret_cast<const __m256i*>(input1)));
        const __m256i result2 = _mm256_xor_si256(
            _mm256_loadu_si256(reinterpret_cast<const __m256i*>(output2)),
            _mm256_loadu_si256(reinterpret_cast<const __m256i*>(input2)));
        const __m256i result3 = _mm256_xor_si256(
            _mm256_loadu_si256(reinterpret_cast<const __m256i*>(output3)),
            _mm256_loadu_si256(reinterpret_cast<const __m256i*>(input3)));
        _mm256_storeu_si256(reinterpret_cast<__m256i*>(output0), result0);
        _mm256_storeu_si256(reinterpret_cast<__m256i*>(output1), result1);
        _mm256_storeu_si256(reinterpret_cast<__m256i*>(output2), result2);
        _mm256_storeu_si256(reinterpret_cast<__m256i*>(output3), result3);
        output0 += 32;
        output1 += 32;
        output2 += 32;
        output3 += 32;
        input0 += 32;
        input1 += 32;
        input2 += 32;
        input3 += 32;
        byte_count -= 32;
    }
    while (byte_count-- != 0)
    {
        *output0++ ^= *input0++;
        *output1++ ^= *input1++;
        *output2++ ^= *input2++;
        *output3++ ^= *input3++;
    }
}

#if defined(_MSC_VER)
#define LEO2_AVX2_FORCE_INLINE __forceinline
#else
#define LEO2_AVX2_FORCE_INLINE inline __attribute__((always_inline))
#endif

#ifdef LEO_HAS_FF8
static void AVX2FF8IFFTButterfly4Nonzero(
    void* value0_pointer, void* value1_pointer,
    void* value2_pointer, void* value3_pointer,
    uint16_t log01, uint16_t log23, uint16_t log02,
    uint64_t byte_count)
{
    uint8_t* value0 = static_cast<uint8_t*>(value0_pointer);
    uint8_t* value1 = static_cast<uint8_t*>(value1_pointer);
    uint8_t* value2 = static_cast<uint8_t*>(value2_pointer);
    uint8_t* value3 = static_cast<uint8_t*>(value3_pointer);
    const __m256i low01 = BroadcastTable(FF8Tables[log01].low);
    const __m256i high01 = BroadcastTable(FF8Tables[log01].high);
    const __m256i low23 = BroadcastTable(FF8Tables[log23].low);
    const __m256i high23 = BroadcastTable(FF8Tables[log23].high);
    const __m256i low02 = BroadcastTable(FF8Tables[log02].low);
    const __m256i high02 = BroadcastTable(FF8Tables[log02].high);

    while (byte_count >= 32)
    {
        __m256i x0 = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(value0));
        __m256i x1 = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(value1));
        __m256i x2 = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(value2));
        __m256i x3 = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(value3));
        x1 = _mm256_xor_si256(x1, x0);
        x0 = _mm256_xor_si256(x0,
            AVX2FF8ProductVector(x1, low01, high01));
        x3 = _mm256_xor_si256(x3, x2);
        x2 = _mm256_xor_si256(x2,
            AVX2FF8ProductVector(x3, low23, high23));
        x2 = _mm256_xor_si256(x2, x0);
        x3 = _mm256_xor_si256(x3, x1);
        x0 = _mm256_xor_si256(x0,
            AVX2FF8ProductVector(x2, low02, high02));
        x1 = _mm256_xor_si256(x1,
            AVX2FF8ProductVector(x3, low02, high02));
        _mm256_storeu_si256(reinterpret_cast<__m256i*>(value0), x0);
        _mm256_storeu_si256(reinterpret_cast<__m256i*>(value1), x1);
        _mm256_storeu_si256(reinterpret_cast<__m256i*>(value2), x2);
        _mm256_storeu_si256(reinterpret_cast<__m256i*>(value3), x3);
        value0 += 32;
        value1 += 32;
        value2 += 32;
        value3 += 32;
        byte_count -= 32;
    }
    while (byte_count-- != 0)
    {
        uint8_t x0 = *value0;
        uint8_t x1 = *value1;
        uint8_t x2 = *value2;
        uint8_t x3 = *value3;
        x1 ^= x0;
        x0 ^= FF8Product(log01, x1);
        x3 ^= x2;
        x2 ^= FF8Product(log23, x3);
        x2 ^= x0;
        x3 ^= x1;
        x0 ^= FF8Product(log02, x2);
        x1 ^= FF8Product(log02, x3);
        *value0++ = x0;
        *value1++ = x1;
        *value2++ = x2;
        *value3++ = x3;
    }
}

static void AVX2FF8IFFTButterfly4Kernel(
    void* value0_pointer, void* value1_pointer,
    void* value2_pointer, void* value3_pointer,
    uint16_t log01, uint16_t log23, uint16_t log02,
    uint64_t byte_count)
{
    static const uint16_t kZeroSkew = 255;
    if (log01 != kZeroSkew && log23 != kZeroSkew && log02 != kZeroSkew)
    {
        AVX2FF8IFFTButterfly4Nonzero(
            value0_pointer, value1_pointer, value2_pointer, value3_pointer,
            log01, log23, log02, byte_count);
        return;
    }

    uint8_t* value0 = static_cast<uint8_t*>(value0_pointer);
    uint8_t* value1 = static_cast<uint8_t*>(value1_pointer);
    uint8_t* value2 = static_cast<uint8_t*>(value2_pointer);
    uint8_t* value3 = static_cast<uint8_t*>(value3_pointer);

    __m256i low01 = _mm256_setzero_si256();
    __m256i high01 = _mm256_setzero_si256();
    __m256i low23 = _mm256_setzero_si256();
    __m256i high23 = _mm256_setzero_si256();
    __m256i low02 = _mm256_setzero_si256();
    __m256i high02 = _mm256_setzero_si256();
    if (log01 != kZeroSkew)
    {
        low01 = BroadcastTable(FF8Tables[log01].low);
        high01 = BroadcastTable(FF8Tables[log01].high);
    }
    if (log23 != kZeroSkew)
    {
        low23 = BroadcastTable(FF8Tables[log23].low);
        high23 = BroadcastTable(FF8Tables[log23].high);
    }
    if (log02 != kZeroSkew)
    {
        low02 = BroadcastTable(FF8Tables[log02].low);
        high02 = BroadcastTable(FF8Tables[log02].high);
    }

    while (byte_count >= 32)
    {
        __m256i x0 = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(value0));
        __m256i x1 = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(value1));
        __m256i x2 = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(value2));
        __m256i x3 = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(value3));
        x1 = _mm256_xor_si256(x1, x0);
        if (log01 != kZeroSkew)
            x0 = _mm256_xor_si256(x0,
                AVX2FF8ProductVector(x1, low01, high01));
        x3 = _mm256_xor_si256(x3, x2);
        if (log23 != kZeroSkew)
            x2 = _mm256_xor_si256(x2,
                AVX2FF8ProductVector(x3, low23, high23));
        x2 = _mm256_xor_si256(x2, x0);
        x3 = _mm256_xor_si256(x3, x1);
        if (log02 != kZeroSkew)
        {
            x0 = _mm256_xor_si256(x0,
                AVX2FF8ProductVector(x2, low02, high02));
            x1 = _mm256_xor_si256(x1,
                AVX2FF8ProductVector(x3, low02, high02));
        }
        _mm256_storeu_si256(reinterpret_cast<__m256i*>(value0), x0);
        _mm256_storeu_si256(reinterpret_cast<__m256i*>(value1), x1);
        _mm256_storeu_si256(reinterpret_cast<__m256i*>(value2), x2);
        _mm256_storeu_si256(reinterpret_cast<__m256i*>(value3), x3);
        value0 += 32;
        value1 += 32;
        value2 += 32;
        value3 += 32;
        byte_count -= 32;
    }
    while (byte_count-- != 0)
    {
        uint8_t x0 = *value0;
        uint8_t x1 = *value1;
        uint8_t x2 = *value2;
        uint8_t x3 = *value3;
        x1 ^= x0;
        if (log01 != kZeroSkew)
            x0 ^= FF8Product(log01, x1);
        x3 ^= x2;
        if (log23 != kZeroSkew)
            x2 ^= FF8Product(log23, x3);
        x2 ^= x0;
        x3 ^= x1;
        if (log02 != kZeroSkew)
        {
            x0 ^= FF8Product(log02, x2);
            x1 ^= FF8Product(log02, x3);
        }
        *value0++ = x0;
        *value1++ = x1;
        *value2++ = x2;
        *value3++ = x3;
    }
}

// The final radix-four layer for every message block after the first feeds an
// XOR accumulator and never observes the temporary work again.  Keep the
// transformed values in registers and accumulate them directly, matching the
// legacy encoder's single memory pass instead of storing work and rereading it
// through AVX2XorMemory4.
template<bool AllNonzero>
static void AVX2FF8IFFTButterfly4XorKernel(
    const void* value0_pointer, const void* value1_pointer,
    const void* value2_pointer, const void* value3_pointer,
    void* output0_pointer, void* output1_pointer,
    void* output2_pointer, void* output3_pointer,
    uint16_t log01, uint16_t log23, uint16_t log02,
    uint64_t byte_count)
{
    static const uint16_t kZeroSkew = 255;
    const uint8_t* value0 = static_cast<const uint8_t*>(value0_pointer);
    const uint8_t* value1 = static_cast<const uint8_t*>(value1_pointer);
    const uint8_t* value2 = static_cast<const uint8_t*>(value2_pointer);
    const uint8_t* value3 = static_cast<const uint8_t*>(value3_pointer);
    uint8_t* output0 = static_cast<uint8_t*>(output0_pointer);
    uint8_t* output1 = static_cast<uint8_t*>(output1_pointer);
    uint8_t* output2 = static_cast<uint8_t*>(output2_pointer);
    uint8_t* output3 = static_cast<uint8_t*>(output3_pointer);

    __m256i low01 = _mm256_setzero_si256();
    __m256i high01 = _mm256_setzero_si256();
    __m256i low23 = _mm256_setzero_si256();
    __m256i high23 = _mm256_setzero_si256();
    __m256i low02 = _mm256_setzero_si256();
    __m256i high02 = _mm256_setzero_si256();
    if (AllNonzero || log01 != kZeroSkew)
    {
        low01 = BroadcastTable(FF8Tables[log01].low);
        high01 = BroadcastTable(FF8Tables[log01].high);
    }
    if (AllNonzero || log23 != kZeroSkew)
    {
        low23 = BroadcastTable(FF8Tables[log23].low);
        high23 = BroadcastTable(FF8Tables[log23].high);
    }
    if (AllNonzero || log02 != kZeroSkew)
    {
        low02 = BroadcastTable(FF8Tables[log02].low);
        high02 = BroadcastTable(FF8Tables[log02].high);
    }

    while (byte_count >= 32)
    {
        __m256i x0 = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(value0));
        __m256i x1 = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(value1));
        __m256i x2 = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(value2));
        __m256i x3 = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(value3));
        x1 = _mm256_xor_si256(x1, x0);
        if (AllNonzero || log01 != kZeroSkew)
            x0 = _mm256_xor_si256(x0,
                AVX2FF8ProductVector(x1, low01, high01));
        x3 = _mm256_xor_si256(x3, x2);
        if (AllNonzero || log23 != kZeroSkew)
            x2 = _mm256_xor_si256(x2,
                AVX2FF8ProductVector(x3, low23, high23));
        x2 = _mm256_xor_si256(x2, x0);
        x3 = _mm256_xor_si256(x3, x1);
        if (AllNonzero || log02 != kZeroSkew)
        {
            x0 = _mm256_xor_si256(x0,
                AVX2FF8ProductVector(x2, low02, high02));
            x1 = _mm256_xor_si256(x1,
                AVX2FF8ProductVector(x3, low02, high02));
        }
        x0 = _mm256_xor_si256(x0, _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(output0)));
        x1 = _mm256_xor_si256(x1, _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(output1)));
        x2 = _mm256_xor_si256(x2, _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(output2)));
        x3 = _mm256_xor_si256(x3, _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(output3)));
        _mm256_storeu_si256(reinterpret_cast<__m256i*>(output0), x0);
        _mm256_storeu_si256(reinterpret_cast<__m256i*>(output1), x1);
        _mm256_storeu_si256(reinterpret_cast<__m256i*>(output2), x2);
        _mm256_storeu_si256(reinterpret_cast<__m256i*>(output3), x3);
        value0 += 32;
        value1 += 32;
        value2 += 32;
        value3 += 32;
        output0 += 32;
        output1 += 32;
        output2 += 32;
        output3 += 32;
        byte_count -= 32;
    }
    while (byte_count-- != 0)
    {
        uint8_t x0 = *value0++;
        uint8_t x1 = *value1++;
        uint8_t x2 = *value2++;
        uint8_t x3 = *value3++;
        x1 ^= x0;
        if (AllNonzero || log01 != kZeroSkew)
            x0 ^= FF8Product(log01, x1);
        x3 ^= x2;
        if (AllNonzero || log23 != kZeroSkew)
            x2 ^= FF8Product(log23, x3);
        x2 ^= x0;
        x3 ^= x1;
        if (AllNonzero || log02 != kZeroSkew)
        {
            x0 ^= FF8Product(log02, x2);
            x1 ^= FF8Product(log02, x3);
        }
        *output0++ ^= x0;
        *output1++ ^= x1;
        *output2++ ^= x2;
        *output3++ ^= x3;
    }
}

static void AVX2FF8IFFTButterfly4(
    void* value0, void* value1, void* value2, void* value3,
    uint16_t log01, uint16_t log23, uint16_t log02,
    uint64_t byte_count)
{
    AVX2FF8IFFTButterfly4Kernel(value0, value1, value2, value3,
        log01, log23, log02, byte_count);
}

static void AVX2FF8FFTButterfly4Fused(
    void* value0_pointer, void* value1_pointer,
    void* value2_pointer, void* value3_pointer,
    uint16_t log01, uint16_t log23, uint16_t log02,
    uint64_t byte_count)
{
    static const uint16_t kZeroSkew = 255;
    uint8_t* value0 = static_cast<uint8_t*>(value0_pointer);
    uint8_t* value1 = static_cast<uint8_t*>(value1_pointer);
    uint8_t* value2 = static_cast<uint8_t*>(value2_pointer);
    uint8_t* value3 = static_cast<uint8_t*>(value3_pointer);
    __m256i low01 = _mm256_setzero_si256();
    __m256i high01 = _mm256_setzero_si256();
    __m256i low23 = _mm256_setzero_si256();
    __m256i high23 = _mm256_setzero_si256();
    __m256i low02 = _mm256_setzero_si256();
    __m256i high02 = _mm256_setzero_si256();
    if (log01 != kZeroSkew)
    {
        low01 = BroadcastTable(FF8Tables[log01].low);
        high01 = BroadcastTable(FF8Tables[log01].high);
    }
    if (log23 != kZeroSkew)
    {
        low23 = BroadcastTable(FF8Tables[log23].low);
        high23 = BroadcastTable(FF8Tables[log23].high);
    }
    if (log02 != kZeroSkew)
    {
        low02 = BroadcastTable(FF8Tables[log02].low);
        high02 = BroadcastTable(FF8Tables[log02].high);
    }

    while (byte_count >= 32)
    {
        __m256i x0 = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(value0));
        __m256i x1 = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(value1));
        __m256i x2 = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(value2));
        __m256i x3 = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(value3));
        if (log02 != kZeroSkew)
        {
            x0 = _mm256_xor_si256(x0,
                AVX2FF8ProductVector(x2, low02, high02));
            x1 = _mm256_xor_si256(x1,
                AVX2FF8ProductVector(x3, low02, high02));
        }
        x2 = _mm256_xor_si256(x2, x0);
        x3 = _mm256_xor_si256(x3, x1);
        if (log01 != kZeroSkew)
            x0 = _mm256_xor_si256(x0,
                AVX2FF8ProductVector(x1, low01, high01));
        x1 = _mm256_xor_si256(x1, x0);
        if (log23 != kZeroSkew)
            x2 = _mm256_xor_si256(x2,
                AVX2FF8ProductVector(x3, low23, high23));
        x3 = _mm256_xor_si256(x3, x2);
        _mm256_storeu_si256(reinterpret_cast<__m256i*>(value0), x0);
        _mm256_storeu_si256(reinterpret_cast<__m256i*>(value1), x1);
        _mm256_storeu_si256(reinterpret_cast<__m256i*>(value2), x2);
        _mm256_storeu_si256(reinterpret_cast<__m256i*>(value3), x3);
        value0 += 32;
        value1 += 32;
        value2 += 32;
        value3 += 32;
        byte_count -= 32;
    }
    while (byte_count-- != 0)
    {
        uint8_t x0 = *value0;
        uint8_t x1 = *value1;
        uint8_t x2 = *value2;
        uint8_t x3 = *value3;
        if (log02 != kZeroSkew)
        {
            x0 ^= FF8Product(log02, x2);
            x1 ^= FF8Product(log02, x3);
        }
        x2 ^= x0;
        x3 ^= x1;
        if (log01 != kZeroSkew)
            x0 ^= FF8Product(log01, x1);
        x1 ^= x0;
        if (log23 != kZeroSkew)
            x2 ^= FF8Product(log23, x3);
        x3 ^= x2;
        *value0++ = x0;
        *value1++ = x1;
        *value2++ = x2;
        *value3++ = x3;
    }
}

static void AVX2FF8FFTButterfly4(
    void* value0, void* value1, void* value2, void* value3,
    uint16_t log01, uint16_t log23, uint16_t log02,
    uint64_t byte_count)
{
    static const uint16_t kZeroSkew = 255;
    static const uint64_t kFusedByteLimit = 1024;
    if (byte_count <= kFusedByteLimit)
    {
        AVX2FF8FFTButterfly4Fused(
            value0, value1, value2, value3,
            log01, log23, log02, byte_count);
        return;
    }
    if (log02 == kZeroSkew)
    {
        AVX2XorMemory(value2, value0, byte_count);
        AVX2XorMemory(value3, value1, byte_count);
    }
    else
    {
        AVX2FF8Butterfly2<false>(value0, value2, log02, byte_count);
        AVX2FF8Butterfly2<false>(value1, value3, log02, byte_count);
    }
    if (log01 == kZeroSkew)
        AVX2XorMemory(value1, value0, byte_count);
    else
        AVX2FF8Butterfly2<false>(value0, value1, log01, byte_count);
    if (log23 == kZeroSkew)
        AVX2XorMemory(value3, value2, byte_count);
    else
        AVX2FF8Butterfly2<false>(value2, value3, log23, byte_count);
}

template<bool Inverse, bool Input3Zero = false>
static void AVX2FF8Butterfly4Out(
    const void* input0_pointer, const void* input1_pointer,
    const void* input2_pointer, const void* input3_pointer,
    void* output0_pointer, void* output1_pointer,
    void* output2_pointer, void* output3_pointer,
    uint16_t log01, uint16_t log23, uint16_t log02,
    uint64_t byte_count)
{
    static const uint16_t kZeroSkew = 255;
    const uint8_t* input0 = static_cast<const uint8_t*>(input0_pointer);
    const uint8_t* input1 = static_cast<const uint8_t*>(input1_pointer);
    const uint8_t* input2 = static_cast<const uint8_t*>(input2_pointer);
    const uint8_t* input3 = static_cast<const uint8_t*>(input3_pointer);
    uint8_t* output0 = static_cast<uint8_t*>(output0_pointer);
    uint8_t* output1 = static_cast<uint8_t*>(output1_pointer);
    uint8_t* output2 = static_cast<uint8_t*>(output2_pointer);
    uint8_t* output3 = static_cast<uint8_t*>(output3_pointer);
    __m256i low01 = _mm256_setzero_si256();
    __m256i high01 = _mm256_setzero_si256();
    __m256i low23 = _mm256_setzero_si256();
    __m256i high23 = _mm256_setzero_si256();
    __m256i low02 = _mm256_setzero_si256();
    __m256i high02 = _mm256_setzero_si256();
    if (log01 != kZeroSkew)
    {
        low01 = BroadcastTable(FF8Tables[log01].low);
        high01 = BroadcastTable(FF8Tables[log01].high);
    }
    if (log23 != kZeroSkew)
    {
        low23 = BroadcastTable(FF8Tables[log23].low);
        high23 = BroadcastTable(FF8Tables[log23].high);
    }
    if (log02 != kZeroSkew)
    {
        low02 = BroadcastTable(FF8Tables[log02].low);
        high02 = BroadcastTable(FF8Tables[log02].high);
    }
    while (byte_count >= 32)
    {
        __m256i x0 = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(input0));
        __m256i x1 = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(input1));
        __m256i x2 = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(input2));
        __m256i x3 = Input3Zero
            ? _mm256_setzero_si256()
            : _mm256_loadu_si256(
                reinterpret_cast<const __m256i*>(input3));
        if (Inverse)
        {
            x1 = _mm256_xor_si256(x1, x0);
            if (log01 != kZeroSkew)
                x0 = _mm256_xor_si256(x0,
                    AVX2FF8ProductVector(x1, low01, high01));
            x3 = _mm256_xor_si256(x3, x2);
            if (log23 != kZeroSkew)
                x2 = _mm256_xor_si256(x2,
                    AVX2FF8ProductVector(x3, low23, high23));
            x2 = _mm256_xor_si256(x2, x0);
            x3 = _mm256_xor_si256(x3, x1);
            if (log02 != kZeroSkew)
            {
                x0 = _mm256_xor_si256(x0,
                    AVX2FF8ProductVector(x2, low02, high02));
                x1 = _mm256_xor_si256(x1,
                    AVX2FF8ProductVector(x3, low02, high02));
            }
        }
        else
        {
            if (log02 != kZeroSkew)
            {
                x0 = _mm256_xor_si256(x0,
                    AVX2FF8ProductVector(x2, low02, high02));
                x1 = _mm256_xor_si256(x1,
                    AVX2FF8ProductVector(x3, low02, high02));
            }
            x2 = _mm256_xor_si256(x2, x0);
            x3 = _mm256_xor_si256(x3, x1);
            if (log01 != kZeroSkew)
                x0 = _mm256_xor_si256(x0,
                    AVX2FF8ProductVector(x1, low01, high01));
            x1 = _mm256_xor_si256(x1, x0);
            if (log23 != kZeroSkew)
                x2 = _mm256_xor_si256(x2,
                    AVX2FF8ProductVector(x3, low23, high23));
            x3 = _mm256_xor_si256(x3, x2);
        }
        _mm256_storeu_si256(reinterpret_cast<__m256i*>(output0), x0);
        _mm256_storeu_si256(reinterpret_cast<__m256i*>(output1), x1);
        _mm256_storeu_si256(reinterpret_cast<__m256i*>(output2), x2);
        _mm256_storeu_si256(reinterpret_cast<__m256i*>(output3), x3);
        input0 += 32;
        input1 += 32;
        input2 += 32;
        if (!Input3Zero)
            input3 += 32;
        output0 += 32;
        output1 += 32;
        output2 += 32;
        output3 += 32;
        byte_count -= 32;
    }
    while (byte_count-- != 0)
    {
        uint8_t x0 = *input0++;
        uint8_t x1 = *input1++;
        uint8_t x2 = *input2++;
        uint8_t x3 = Input3Zero ? 0 : *input3++;
        if (Inverse)
        {
            x1 ^= x0;
            if (log01 != kZeroSkew)
                x0 ^= FF8Product(log01, x1);
            x3 ^= x2;
            if (log23 != kZeroSkew)
                x2 ^= FF8Product(log23, x3);
            x2 ^= x0;
            x3 ^= x1;
            if (log02 != kZeroSkew)
            {
                x0 ^= FF8Product(log02, x2);
                x1 ^= FF8Product(log02, x3);
            }
        }
        else
        {
            if (log02 != kZeroSkew)
            {
                x0 ^= FF8Product(log02, x2);
                x1 ^= FF8Product(log02, x3);
            }
            x2 ^= x0;
            x3 ^= x1;
            if (log01 != kZeroSkew)
                x0 ^= FF8Product(log01, x1);
            x1 ^= x0;
            if (log23 != kZeroSkew)
                x2 ^= FF8Product(log23, x3);
            x3 ^= x2;
        }
        *output0++ = x0;
        *output1++ = x1;
        *output2++ = x2;
        *output3++ = x3;
    }
}

static void AVX2FF8IFFTButterfly4Out(
    const void* input0, const void* input1,
    const void* input2, const void* input3,
    void* output0, void* output1, void* output2, void* output3,
    uint16_t log01, uint16_t log23, uint16_t log02,
    uint64_t byte_count)
{
    AVX2FF8Butterfly4Out<true>(
        input0, input1, input2, input3,
        output0, output1, output2, output3,
        log01, log23, log02, byte_count);
}

#if defined(LEO2_AVX512_VARIANT)

#if defined(_MSC_VER)
# define LEO2_AVX2_NOINLINE __declspec(noinline)
#elif defined(__GNUC__) || defined(__clang__)
# define LEO2_AVX2_NOINLINE __attribute__((noinline))
#else
# define LEO2_AVX2_NOINLINE
#endif

static LEO2_AVX2_NOINLINE void AVX2FF8IFFTButterfly4OutLastZero(
    const void* input0, const void* input1, const void* input2,
    void* output0, void* output1, void* output2, void* output3,
    uint16_t log01, uint16_t log23, uint16_t log02,
    uint64_t byte_count)
{
    AVX2FF8Butterfly4Out<true, true>(
        input0, input1, input2, NULL,
        output0, output1, output2, output3,
        log01, log23, log02, byte_count);
}

#undef LEO2_AVX2_NOINLINE

#endif // LEO2_AVX512_VARIANT

static void AVX2FF8WeightedIFFTButterfly4(
    const void* input0_pointer, const void* input1_pointer,
    const void* input2_pointer, const void* input3_pointer,
    void* output0_pointer, void* output1_pointer,
    void* output2_pointer, void* output3_pointer,
    uint16_t weight_log0, uint16_t weight_log1,
    uint16_t weight_log2, uint16_t weight_log3,
    uint8_t live_mask,
    uint16_t log01, uint16_t log23, uint16_t log02,
    uint64_t byte_count)
{
    static const uint16_t kModulus = 255;
    const uint8_t* inputs[4] = {
        static_cast<const uint8_t*>(input0_pointer),
        static_cast<const uint8_t*>(input1_pointer),
        static_cast<const uint8_t*>(input2_pointer),
        static_cast<const uint8_t*>(input3_pointer)
    };
    uint8_t* outputs[4] = {
        static_cast<uint8_t*>(output0_pointer),
        static_cast<uint8_t*>(output1_pointer),
        static_cast<uint8_t*>(output2_pointer),
        static_cast<uint8_t*>(output3_pointer)
    };
#ifdef LEO_DEBUG
    const void* alias_inputs[4] = {
        input0_pointer, input1_pointer, input2_pointer, input3_pointer
    };
    const void* alias_outputs[4] = {
        output0_pointer, output1_pointer, output2_pointer, output3_pointer
    };
    LEO_DEBUG_ASSERT(IsWeightedIFFTButterfly4AliasingValid(
        alias_inputs, alias_outputs, live_mask, byte_count));
#endif
    const uint16_t weights[4] = {
        weight_log0, weight_log1, weight_log2, weight_log3
    };
    __m256i low01 = _mm256_setzero_si256();
    __m256i high01 = _mm256_setzero_si256();
    __m256i low23 = _mm256_setzero_si256();
    __m256i high23 = _mm256_setzero_si256();
    __m256i low02 = _mm256_setzero_si256();
    __m256i high02 = _mm256_setzero_si256();
    if (log01 != kModulus)
    {
        low01 = BroadcastTable(FF8Tables[log01].low);
        high01 = BroadcastTable(FF8Tables[log01].high);
    }
    if (log23 != kModulus)
    {
        low23 = BroadcastTable(FF8Tables[log23].low);
        high23 = BroadcastTable(FF8Tables[log23].high);
    }
    if (log02 != kModulus)
    {
        low02 = BroadcastTable(FF8Tables[log02].low);
        high02 = BroadcastTable(FF8Tables[log02].high);
    }
    uint64_t offset = 0;
    while (byte_count - offset >= 32)
    {
        __m256i x[4];
        for (unsigned lane = 0; lane < 4; ++lane)
        {
            if ((live_mask & (1U << lane)) != 0)
            {
                x[lane] = _mm256_loadu_si256(
                    reinterpret_cast<const __m256i*>(inputs[lane] + offset));
                x[lane] = AVX2FF8ApplyWeight(x[lane], weights[lane]);
            }
            else
                x[lane] = _mm256_setzero_si256();
        }
        x[1] = _mm256_xor_si256(x[1], x[0]);
        if (log01 != kModulus)
            x[0] = _mm256_xor_si256(x[0],
                AVX2FF8ProductVector(x[1], low01, high01));
        x[3] = _mm256_xor_si256(x[3], x[2]);
        if (log23 != kModulus)
            x[2] = _mm256_xor_si256(x[2],
                AVX2FF8ProductVector(x[3], low23, high23));
        x[2] = _mm256_xor_si256(x[2], x[0]);
        x[3] = _mm256_xor_si256(x[3], x[1]);
        if (log02 != kModulus)
        {
            x[0] = _mm256_xor_si256(x[0],
                AVX2FF8ProductVector(x[2], low02, high02));
            x[1] = _mm256_xor_si256(x[1],
                AVX2FF8ProductVector(x[3], low02, high02));
        }
        for (unsigned lane = 0; lane < 4; ++lane)
            _mm256_storeu_si256(
                reinterpret_cast<__m256i*>(outputs[lane] + offset), x[lane]);
        offset += 32;
    }
    while (offset < byte_count)
    {
        uint8_t x[4] = {};
        for (unsigned lane = 0; lane < 4; ++lane)
        {
            if ((live_mask & (1U << lane)) == 0)
                continue;
            x[lane] = inputs[lane][offset];
            if (weights[lane] != 0 && weights[lane] != kModulus)
                x[lane] = FF8Product(weights[lane], x[lane]);
        }
        x[1] ^= x[0];
        if (log01 != kModulus)
            x[0] ^= FF8Product(log01, x[1]);
        x[3] ^= x[2];
        if (log23 != kModulus)
            x[2] ^= FF8Product(log23, x[3]);
        x[2] ^= x[0];
        x[3] ^= x[1];
        if (log02 != kModulus)
        {
            x[0] ^= FF8Product(log02, x[2]);
            x[1] ^= FF8Product(log02, x[3]);
        }
        for (unsigned lane = 0; lane < 4; ++lane)
            outputs[lane][offset] = x[lane];
        ++offset;
    }
}

static void AVX2FF8FFTButterfly4Out(
    const void* input0, const void* input1,
    const void* input2, const void* input3,
    void* output0, void* output1, void* output2, void* output3,
    uint16_t log01, uint16_t log23, uint16_t log02,
    uint64_t byte_count)
{
    AVX2FF8Butterfly4Out<false>(
        input0, input1, input2, input3,
        output0, output1, output2, output3,
        log01, log23, log02, byte_count);
}

static void AVX2FF8IFFTButterfly4Range(
    void* const* work, unsigned distance,
    uint16_t log01, uint16_t log23, uint16_t log02,
    uint64_t byte_count, bool prefer_fused)
{
    (void)prefer_fused;
    for (unsigned i = 0; i < distance; ++i)
    {
        AVX2FF8IFFTButterfly4Kernel(
            work[i], work[i + distance],
            work[i + distance * 2U], work[i + distance * 3U],
            log01, log23, log02, byte_count);
    }
}

static void AVX2FF8FFTButterfly4Range(
    void* const* work, unsigned distance,
    uint16_t log01, uint16_t log23, uint16_t log02,
    uint64_t byte_count, bool prefer_fused)
{
    if (prefer_fused)
    {
        for (unsigned i = 0; i < distance; ++i)
        {
            AVX2FF8FFTButterfly4Fused(
                work[i], work[i + distance],
                work[i + distance * 2U], work[i + distance * 3U],
                log01, log23, log02, byte_count);
        }
        return;
    }
    for (unsigned i = 0; i < distance; ++i)
    {
        AVX2FF8FFTButterfly4(
            work[i], work[i + distance],
            work[i + distance * 2U], work[i + distance * 3U],
            log01, log23, log02, byte_count);
    }
}

static void AVX2FF8IFFTButterfly4XorRange(
    void* const* work, void* const* xor_output, unsigned distance,
    uint16_t log01, uint16_t log23, uint16_t log02,
    uint64_t byte_count)
{
    static const uint16_t kZeroSkew = 255;
    const bool all_nonzero =
        log01 != kZeroSkew && log23 != kZeroSkew && log02 != kZeroSkew;
    if (all_nonzero)
    {
        for (unsigned i = 0; i < distance; ++i)
            AVX2FF8IFFTButterfly4XorKernel<true>(
                work[i], work[i + distance],
                work[i + distance * 2U], work[i + distance * 3U],
                xor_output[i], xor_output[i + distance],
                xor_output[i + distance * 2U],
                xor_output[i + distance * 3U],
                log01, log23, log02, byte_count);
        return;
    }
    for (unsigned i = 0; i < distance; ++i)
        AVX2FF8IFFTButterfly4XorKernel<false>(
            work[i], work[i + distance],
            work[i + distance * 2U], work[i + distance * 3U],
            xor_output[i], xor_output[i + distance],
            xor_output[i + distance * 2U],
            xor_output[i + distance * 3U],
            log01, log23, log02, byte_count);
}

#if !defined(LEO2_AVX512_VARIANT)
static void AVX2FF8HighEncodeSmall(
    const void* const* data,
    uint32_t original_count,
    void* const* work,
    uint32_t side,
    const uint8_t* inverse_skew,
    const uint8_t* forward_skew,
    uint64_t byte_count)
{
    static const uint16_t kZeroSkew = 255;

    // Materialize the first complete message block directly into the
    // accumulator.  This avoids both a separate output clear and the
    // redundant zero-valued accumulator load in the first XOR transform.
    if (side == 2)
    {
        AVX2FF8IFFTButterfly2Out(
            data[0], data[1], work[0], work[1],
            inverse_skew[1], byte_count);
    }
    else
    {
        AVX2FF8IFFTButterfly4Out(
            data[0], data[1], data[2], data[3],
            work[0], work[1], work[2], work[3],
            inverse_skew[1], inverse_skew[3], inverse_skew[2],
            byte_count);
    }

    for (uint32_t base = side; base < original_count; base += side)
    {
        const uint32_t remaining = std::min(
            side, original_count - base);
        const void* source[4] = { NULL, NULL, NULL, NULL };
        if (remaining == side)
        {
            for (uint32_t lane = 0; lane < side; ++lane)
                source[lane] = data[base + lane];
        }
        else
        {
            // work[side..2*side) is the established encoder temporary half.
            // A high-rate block can have at most one shortened suffix, so the
            // only retained copy/zero pass is bounded independently of K.
            for (uint32_t lane = 0; lane < side; ++lane)
            {
                if (lane < remaining)
                    std::memcpy(work[side + lane], data[base + lane],
                        static_cast<size_t>(byte_count));
                else
                    std::memset(work[side + lane], 0,
                        static_cast<size_t>(byte_count));
                source[lane] = work[side + lane];
            }
        }

        if (side == 2)
        {
            const uint16_t log = inverse_skew[base + 1U];
            if (log == kZeroSkew)
            {
                AVX2XorMemory(work[0], source[0], byte_count);
                AVX2XorMemory2To1(
                    work[1], source[0], source[1], byte_count);
            }
            else
            {
                AVX2FF8IFFTButterfly2Xor(
                    source[0], source[1], work[0], work[1],
                    log, byte_count);
            }
        }
        else
        {
            const uint16_t log01 = inverse_skew[base + 1U];
            const uint16_t log23 = inverse_skew[base + 3U];
            const uint16_t log02 = inverse_skew[base + 2U];
            if (log01 != kZeroSkew && log23 != kZeroSkew &&
                log02 != kZeroSkew)
            {
                AVX2FF8IFFTButterfly4XorKernel<true>(
                    source[0], source[1], source[2], source[3],
                    work[0], work[1], work[2], work[3],
                    log01, log23, log02, byte_count);
            }
            else
            {
                AVX2FF8IFFTButterfly4XorKernel<false>(
                    source[0], source[1], source[2], source[3],
                    work[0], work[1], work[2], work[3],
                    log01, log23, log02, byte_count);
            }
        }
    }

    if (side == 2)
    {
        const uint16_t log = forward_skew[1];
        if (log == kZeroSkew)
            AVX2XorMemory(work[1], work[0], byte_count);
        else
            AVX2FF8FFTButterfly2(work[0], work[1], log, byte_count);
    }
    else
    {
        AVX2FF8FFTButterfly4(
            work[0], work[1], work[2], work[3],
            forward_skew[1], forward_skew[3], forward_skew[2],
            byte_count);
    }
}
#endif

#if defined(LEO2_AVX512_VARIANT)

struct AVX2FF8T8Pair
{
    __m256i low;
    __m256i high;
};

static LEO2_AVX2_FORCE_INLINE void AVX2FF8T8Xor(
    AVX2FF8T8Pair& destination,
    const AVX2FF8T8Pair& source)
{
    destination.low = _mm256_xor_si256(destination.low, source.low);
    destination.high = _mm256_xor_si256(destination.high, source.high);
}

static LEO2_AVX2_FORCE_INLINE void AVX2FF8T8MulAdd(
    AVX2FF8T8Pair& destination,
    const AVX2FF8T8Pair& source,
    uint16_t log)
{
    const FF8NibbleTable& table = FF8Tables[log];
    const __m256i low_table = BroadcastTable(table.low);
    const __m256i high_table = BroadcastTable(table.high);
    destination.low = _mm256_xor_si256(destination.low,
        AVX2FF8ProductVector(source.low, low_table, high_table));
    destination.high = _mm256_xor_si256(destination.high,
        AVX2FF8ProductVector(source.high, low_table, high_table));
}

static LEO2_AVX2_FORCE_INLINE void AVX2FF8T8MulAdd2(
    AVX2FF8T8Pair& destination0,
    const AVX2FF8T8Pair& source0,
    AVX2FF8T8Pair& destination1,
    const AVX2FF8T8Pair& source1,
    uint16_t log)
{
    const FF8NibbleTable& table = FF8Tables[log];
    const __m256i low_table = BroadcastTable(table.low);
    const __m256i high_table = BroadcastTable(table.high);
    destination0.low = _mm256_xor_si256(destination0.low,
        AVX2FF8ProductVector(source0.low, low_table, high_table));
    destination0.high = _mm256_xor_si256(destination0.high,
        AVX2FF8ProductVector(source0.high, low_table, high_table));
    destination1.low = _mm256_xor_si256(destination1.low,
        AVX2FF8ProductVector(source1.low, low_table, high_table));
    destination1.high = _mm256_xor_si256(destination1.high,
        AVX2FF8ProductVector(source1.high, low_table, high_table));
}

static LEO2_AVX2_FORCE_INLINE void AVX2FF8T8MulAdd4(
    AVX2FF8T8Pair& destination0,
    const AVX2FF8T8Pair& source0,
    AVX2FF8T8Pair& destination1,
    const AVX2FF8T8Pair& source1,
    AVX2FF8T8Pair& destination2,
    const AVX2FF8T8Pair& source2,
    AVX2FF8T8Pair& destination3,
    const AVX2FF8T8Pair& source3,
    uint16_t log)
{
    const FF8NibbleTable& table = FF8Tables[log];
    const __m256i low_table = BroadcastTable(table.low);
    const __m256i high_table = BroadcastTable(table.high);
    destination0.low = _mm256_xor_si256(destination0.low,
        AVX2FF8ProductVector(source0.low, low_table, high_table));
    destination0.high = _mm256_xor_si256(destination0.high,
        AVX2FF8ProductVector(source0.high, low_table, high_table));
    destination1.low = _mm256_xor_si256(destination1.low,
        AVX2FF8ProductVector(source1.low, low_table, high_table));
    destination1.high = _mm256_xor_si256(destination1.high,
        AVX2FF8ProductVector(source1.high, low_table, high_table));
    destination2.low = _mm256_xor_si256(destination2.low,
        AVX2FF8ProductVector(source2.low, low_table, high_table));
    destination2.high = _mm256_xor_si256(destination2.high,
        AVX2FF8ProductVector(source2.high, low_table, high_table));
    destination3.low = _mm256_xor_si256(destination3.low,
        AVX2FF8ProductVector(source3.low, low_table, high_table));
    destination3.high = _mm256_xor_si256(destination3.high,
        AVX2FF8ProductVector(source3.high, low_table, high_table));
}

static LEO2_AVX2_FORCE_INLINE void AVX2FF8T8IFFTRadix4(
    AVX2FF8T8Pair& value0,
    AVX2FF8T8Pair& value1,
    AVX2FF8T8Pair& value2,
    AVX2FF8T8Pair& value3,
    uint16_t log01,
    uint16_t log23,
    uint16_t log02)
{
    static const uint16_t kZeroSkew = 255;
    AVX2FF8T8Xor(value1, value0);
    if (log01 != kZeroSkew)
        AVX2FF8T8MulAdd(value0, value1, log01);
    AVX2FF8T8Xor(value3, value2);
    if (log23 != kZeroSkew)
        AVX2FF8T8MulAdd(value2, value3, log23);
    AVX2FF8T8Xor(value2, value0);
    AVX2FF8T8Xor(value3, value1);
    if (log02 != kZeroSkew)
        AVX2FF8T8MulAdd2(value0, value2, value1, value3, log02);
}

static LEO2_AVX2_FORCE_INLINE void AVX2FF8T8IFFTDistance4(
    AVX2FF8T8Pair values[8],
    uint16_t log)
{
    static const uint16_t kZeroSkew = 255;
    AVX2FF8T8Xor(values[4], values[0]);
    AVX2FF8T8Xor(values[5], values[1]);
    AVX2FF8T8Xor(values[6], values[2]);
    AVX2FF8T8Xor(values[7], values[3]);
    if (log != kZeroSkew)
        AVX2FF8T8MulAdd4(
            values[0], values[4], values[1], values[5],
            values[2], values[6], values[3], values[7], log);
}

static LEO2_AVX2_FORCE_INLINE void AVX2FF8T8FFTRadix4Distance2(
    AVX2FF8T8Pair values[8],
    uint16_t log01,
    uint16_t log23,
    uint16_t log02)
{
    static const uint16_t kZeroSkew = 255;
    if (log02 != kZeroSkew)
        AVX2FF8T8MulAdd4(
            values[0], values[4], values[1], values[5],
            values[2], values[6], values[3], values[7], log02);
    AVX2FF8T8Xor(values[4], values[0]);
    AVX2FF8T8Xor(values[5], values[1]);
    AVX2FF8T8Xor(values[6], values[2]);
    AVX2FF8T8Xor(values[7], values[3]);
    if (log01 != kZeroSkew)
        AVX2FF8T8MulAdd2(
            values[0], values[2], values[1], values[3], log01);
    AVX2FF8T8Xor(values[2], values[0]);
    AVX2FF8T8Xor(values[3], values[1]);
    if (log23 != kZeroSkew)
        AVX2FF8T8MulAdd2(
            values[4], values[6], values[5], values[7], log23);
    AVX2FF8T8Xor(values[6], values[4]);
    AVX2FF8T8Xor(values[7], values[5]);
}

static LEO2_AVX2_FORCE_INLINE void AVX2FF8T8FFT2(
    AVX2FF8T8Pair& value0,
    AVX2FF8T8Pair& value1,
    uint16_t log)
{
    static const uint16_t kZeroSkew = 255;
    if (log != kZeroSkew)
        AVX2FF8T8MulAdd(value0, value1, log);
    AVX2FF8T8Xor(value1, value0);
}

static LEO2_AVX2_FORCE_INLINE void AVX2FF8T8ScalarIFFT2(
    uint8_t& value0,
    uint8_t& value1,
    uint16_t log)
{
    value1 ^= value0;
    if (log != 255)
        value0 ^= FF8Product(log, value1);
}

static LEO2_AVX2_FORCE_INLINE void AVX2FF8T8ScalarFFT2(
    uint8_t& value0,
    uint8_t& value1,
    uint16_t log)
{
    if (log != 255)
        value0 ^= FF8Product(log, value1);
    value1 ^= value0;
}

#if defined(_MSC_VER)
# define LEO2_AVX2_T8_NOINLINE __declspec(noinline)
#elif defined(__GNUC__) || defined(__clang__)
# define LEO2_AVX2_T8_NOINLINE __attribute__((noinline))
#else
# define LEO2_AVX2_T8_NOINLINE
#endif

static LEO2_AVX2_T8_NOINLINE void AVX2FF8HighEncodeT8ScalarTail(
    const uint8_t* const input[8],
    uint8_t* const output[8],
    uint64_t lane7_offset_mask,
    const uint8_t* inverse_skew,
    const uint8_t* forward_skew,
    uint64_t offset,
    uint64_t byte_count)
{
    while (offset < byte_count)
    {
        uint8_t values[8];
        for (unsigned lane = 0; lane < 7; ++lane)
            values[lane] = input[lane][offset];
        values[7] = input[7][offset & lane7_offset_mask];

        AVX2FF8T8ScalarIFFT2(values[0], values[1], inverse_skew[1]);
        AVX2FF8T8ScalarIFFT2(values[2], values[3], inverse_skew[3]);
        AVX2FF8T8ScalarIFFT2(values[0], values[2], inverse_skew[2]);
        AVX2FF8T8ScalarIFFT2(values[1], values[3], inverse_skew[2]);
        AVX2FF8T8ScalarIFFT2(values[4], values[5], inverse_skew[5]);
        AVX2FF8T8ScalarIFFT2(values[6], values[7], inverse_skew[7]);
        AVX2FF8T8ScalarIFFT2(values[4], values[6], inverse_skew[6]);
        AVX2FF8T8ScalarIFFT2(values[5], values[7], inverse_skew[6]);
        for (unsigned lane = 0; lane < 4; ++lane)
            AVX2FF8T8ScalarIFFT2(
                values[lane], values[lane + 4], inverse_skew[4]);

        AVX2FF8T8ScalarFFT2(values[0], values[4], forward_skew[4]);
        AVX2FF8T8ScalarFFT2(values[1], values[5], forward_skew[4]);
        AVX2FF8T8ScalarFFT2(values[2], values[6], forward_skew[4]);
        AVX2FF8T8ScalarFFT2(values[3], values[7], forward_skew[4]);
        AVX2FF8T8ScalarFFT2(values[0], values[2], forward_skew[2]);
        AVX2FF8T8ScalarFFT2(values[1], values[3], forward_skew[2]);
        AVX2FF8T8ScalarFFT2(values[4], values[6], forward_skew[6]);
        AVX2FF8T8ScalarFFT2(values[5], values[7], forward_skew[6]);
        AVX2FF8T8ScalarFFT2(values[0], values[1], forward_skew[1]);
        AVX2FF8T8ScalarFFT2(values[2], values[3], forward_skew[3]);
        AVX2FF8T8ScalarFFT2(values[4], values[5], forward_skew[5]);
        AVX2FF8T8ScalarFFT2(values[6], values[7], forward_skew[7]);

        for (unsigned lane = 0; lane < 8; ++lane)
            output[lane][offset] = values[lane];
        ++offset;
    }
}

#undef LEO2_AVX2_T8_NOINLINE

static void AVX2FF8HighEncodeT8(
    const void* const* data,
    void* const* work,
    bool shortened,
    const uint8_t* inverse_skew,
    const uint8_t* forward_skew,
    uint64_t byte_count)
{
    alignas(32) static const uint8_t kZeroShard[64] = {};
    const uint8_t* input[8];
    uint8_t* output[8];
    for (unsigned lane = 0; lane < 7; ++lane)
    {
        input[lane] = static_cast<const uint8_t*>(data[lane]);
        output[lane] = static_cast<uint8_t*>(work[lane]);
    }
    input[7] = shortened
        ? kZeroShard : static_cast<const uint8_t*>(data[7]);
    output[7] = static_cast<uint8_t*>(work[7]);
    const uint64_t lane7_offset_mask = shortened ? 0 : ~uint64_t(0);

    uint64_t offset = 0;
    while (byte_count - offset >= 64)
    {
        AVX2FF8T8Pair values[8];
        for (unsigned lane = 0; lane < 7; ++lane)
        {
            values[lane].low = _mm256_loadu_si256(
                reinterpret_cast<const __m256i*>(input[lane] + offset));
            values[lane].high = _mm256_loadu_si256(
                reinterpret_cast<const __m256i*>(
                    input[lane] + offset + 32));
        }
        const uint64_t lane7_offset = offset & lane7_offset_mask;
        values[7].low = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(input[7] + lane7_offset));
        values[7].high = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(
                input[7] + lane7_offset + 32));

        AVX2FF8T8IFFTRadix4(
            values[0], values[1], values[2], values[3],
            inverse_skew[1], inverse_skew[3], inverse_skew[2]);
        AVX2FF8T8IFFTRadix4(
            values[4], values[5], values[6], values[7],
            inverse_skew[5], inverse_skew[7], inverse_skew[6]);
        AVX2FF8T8IFFTDistance4(values, inverse_skew[4]);

        AVX2FF8T8FFTRadix4Distance2(values,
            forward_skew[2], forward_skew[6], forward_skew[4]);
        AVX2FF8T8FFT2(values[0], values[1], forward_skew[1]);
        AVX2FF8T8FFT2(values[2], values[3], forward_skew[3]);
        AVX2FF8T8FFT2(values[4], values[5], forward_skew[5]);
        AVX2FF8T8FFT2(values[6], values[7], forward_skew[7]);

        for (unsigned lane = 0; lane < 8; ++lane)
        {
            _mm256_storeu_si256(
                reinterpret_cast<__m256i*>(output[lane] + offset),
                values[lane].low);
            _mm256_storeu_si256(
                reinterpret_cast<__m256i*>(output[lane] + offset + 32),
                values[lane].high);
        }
        offset += 64;
    }

    if (offset < byte_count)
        AVX2FF8HighEncodeT8ScalarTail(
            input, output, lane7_offset_mask,
            inverse_skew, forward_skew, offset, byte_count);
}

template<bool Inverse, bool AllNonzero>
static void AVX2FF8Butterfly4RangePreparedImpl(
    void* const* work,
    uint32_t distance,
    uint16_t log01,
    uint16_t log23,
    uint16_t log02,
    uint64_t byte_count)
{
    static const uint16_t kZeroSkew = 255;
    __m256i low01 = _mm256_setzero_si256();
    __m256i high01 = _mm256_setzero_si256();
    __m256i low23 = _mm256_setzero_si256();
    __m256i high23 = _mm256_setzero_si256();
    __m256i low02 = _mm256_setzero_si256();
    __m256i high02 = _mm256_setzero_si256();
    if (AllNonzero || log01 != kZeroSkew)
    {
        low01 = BroadcastTable(FF8Tables[log01].low);
        high01 = BroadcastTable(FF8Tables[log01].high);
    }
    if (AllNonzero || log23 != kZeroSkew)
    {
        low23 = BroadcastTable(FF8Tables[log23].low);
        high23 = BroadcastTable(FF8Tables[log23].high);
    }
    if (AllNonzero || log02 != kZeroSkew)
    {
        low02 = BroadcastTable(FF8Tables[log02].low);
        high02 = BroadcastTable(FF8Tables[log02].high);
    }

    for (uint32_t lane = 0; lane < distance; ++lane)
    {
        uint8_t* value0 = static_cast<uint8_t*>(work[lane]);
        uint8_t* value1 = static_cast<uint8_t*>(work[lane + distance]);
        uint8_t* value2 = static_cast<uint8_t*>(work[lane + distance * 2U]);
        uint8_t* value3 = static_cast<uint8_t*>(work[lane + distance * 3U]);
        uint64_t remaining = byte_count;
        while (remaining >= 32)
        {
            __m256i x0 = _mm256_loadu_si256(
                reinterpret_cast<const __m256i*>(value0));
            __m256i x1 = _mm256_loadu_si256(
                reinterpret_cast<const __m256i*>(value1));
            __m256i x2 = _mm256_loadu_si256(
                reinterpret_cast<const __m256i*>(value2));
            __m256i x3 = _mm256_loadu_si256(
                reinterpret_cast<const __m256i*>(value3));
            if (Inverse)
            {
                x1 = _mm256_xor_si256(x1, x0);
                if (AllNonzero || log01 != kZeroSkew)
                    x0 = _mm256_xor_si256(x0,
                        AVX2FF8ProductVector(x1, low01, high01));
                x3 = _mm256_xor_si256(x3, x2);
                if (AllNonzero || log23 != kZeroSkew)
                    x2 = _mm256_xor_si256(x2,
                        AVX2FF8ProductVector(x3, low23, high23));
                x2 = _mm256_xor_si256(x2, x0);
                x3 = _mm256_xor_si256(x3, x1);
                if (AllNonzero || log02 != kZeroSkew)
                {
                    x0 = _mm256_xor_si256(x0,
                        AVX2FF8ProductVector(x2, low02, high02));
                    x1 = _mm256_xor_si256(x1,
                        AVX2FF8ProductVector(x3, low02, high02));
                }
            }
            else
            {
                if (AllNonzero || log02 != kZeroSkew)
                {
                    x0 = _mm256_xor_si256(x0,
                        AVX2FF8ProductVector(x2, low02, high02));
                    x1 = _mm256_xor_si256(x1,
                        AVX2FF8ProductVector(x3, low02, high02));
                }
                x2 = _mm256_xor_si256(x2, x0);
                x3 = _mm256_xor_si256(x3, x1);
                if (AllNonzero || log01 != kZeroSkew)
                    x0 = _mm256_xor_si256(x0,
                        AVX2FF8ProductVector(x1, low01, high01));
                x1 = _mm256_xor_si256(x1, x0);
                if (AllNonzero || log23 != kZeroSkew)
                    x2 = _mm256_xor_si256(x2,
                        AVX2FF8ProductVector(x3, low23, high23));
                x3 = _mm256_xor_si256(x3, x2);
            }
            _mm256_storeu_si256(reinterpret_cast<__m256i*>(value0), x0);
            _mm256_storeu_si256(reinterpret_cast<__m256i*>(value1), x1);
            _mm256_storeu_si256(reinterpret_cast<__m256i*>(value2), x2);
            _mm256_storeu_si256(reinterpret_cast<__m256i*>(value3), x3);
            value0 += 32;
            value1 += 32;
            value2 += 32;
            value3 += 32;
            remaining -= 32;
        }
        while (remaining-- != 0)
        {
            uint8_t x0 = *value0;
            uint8_t x1 = *value1;
            uint8_t x2 = *value2;
            uint8_t x3 = *value3;
            if (Inverse)
            {
                x1 ^= x0;
                if (AllNonzero || log01 != kZeroSkew)
                    x0 ^= FF8Product(log01, x1);
                x3 ^= x2;
                if (AllNonzero || log23 != kZeroSkew)
                    x2 ^= FF8Product(log23, x3);
                x2 ^= x0;
                x3 ^= x1;
                if (AllNonzero || log02 != kZeroSkew)
                {
                    x0 ^= FF8Product(log02, x2);
                    x1 ^= FF8Product(log02, x3);
                }
            }
            else
            {
                if (AllNonzero || log02 != kZeroSkew)
                {
                    x0 ^= FF8Product(log02, x2);
                    x1 ^= FF8Product(log02, x3);
                }
                x2 ^= x0;
                x3 ^= x1;
                if (AllNonzero || log01 != kZeroSkew)
                    x0 ^= FF8Product(log01, x1);
                x1 ^= x0;
                if (AllNonzero || log23 != kZeroSkew)
                    x2 ^= FF8Product(log23, x3);
                x3 ^= x2;
            }
            *value0++ = x0;
            *value1++ = x1;
            *value2++ = x2;
            *value3++ = x3;
        }
    }
}

template<bool Inverse>
static void AVX2FF8Butterfly4RangePrepared(
    void* const* work,
    uint32_t distance,
    uint16_t log01,
    uint16_t log23,
    uint16_t log02,
    uint64_t byte_count)
{
    static const uint16_t kZeroSkew = 255;
    if (log01 != kZeroSkew && log23 != kZeroSkew && log02 != kZeroSkew)
        AVX2FF8Butterfly4RangePreparedImpl<Inverse, true>(
            work, distance, log01, log23, log02, byte_count);
    else
        AVX2FF8Butterfly4RangePreparedImpl<Inverse, false>(
            work, distance, log01, log23, log02, byte_count);
}

static void AVX2FF8IFFTButterfly2RangePrepared(
    void* const* work,
    uint32_t distance,
    uint16_t log,
    uint64_t byte_count)
{
    static const uint16_t kZeroSkew = 255;
    if (log == kZeroSkew)
    {
        for (uint32_t lane = 0; lane < distance; ++lane)
            AVX2XorMemory(
                work[lane + distance], work[lane], byte_count);
        return;
    }
    const __m256i low = BroadcastTable(FF8Tables[log].low);
    const __m256i high = BroadcastTable(FF8Tables[log].high);
    for (uint32_t lane = 0; lane < distance; ++lane)
    {
        uint8_t* x = static_cast<uint8_t*>(work[lane]);
        uint8_t* y = static_cast<uint8_t*>(work[lane + distance]);
        uint64_t remaining = byte_count;
        while (remaining >= 32)
        {
            __m256i x_value = _mm256_loadu_si256(
                reinterpret_cast<const __m256i*>(x));
            __m256i y_value = _mm256_loadu_si256(
                reinterpret_cast<const __m256i*>(y));
            y_value = _mm256_xor_si256(y_value, x_value);
            x_value = _mm256_xor_si256(x_value,
                AVX2FF8ProductVector(y_value, low, high));
            _mm256_storeu_si256(reinterpret_cast<__m256i*>(x), x_value);
            _mm256_storeu_si256(reinterpret_cast<__m256i*>(y), y_value);
            x += 32;
            y += 32;
            remaining -= 32;
        }
        while (remaining-- != 0)
        {
            *y ^= *x;
            *x ^= FF8Product(log, *y);
            ++x;
            ++y;
        }
    }
}

static void AVX2FF8HighEncodeOneBlock(
    const void* const* data,
    void* const* work,
    uint32_t side_and_flags,
    const uint8_t* inverse_skew,
    const uint8_t* forward_skew,
    uint64_t byte_count)
{
    static const uint16_t kZeroSkew = 255;
    const bool shortened =
        (side_and_flags & kFF8HighEncodeShortenedInput) != 0;
    const uint32_t side =
        side_and_flags & ~kFF8HighEncodeShortenedInput;
    const uint32_t data_count = side - static_cast<uint32_t>(shortened);
    LEO_DEBUG_ASSERT(side == 8 || side == 16 || side == 32 || side == 64);

    if (side == 8)
    {
        AVX2FF8HighEncodeT8(
            data, work, shortened, inverse_skew, forward_skew, byte_count);
        return;
    }

    // Consume the systematic shards directly through the first two inverse
    // layers, exactly as IFFT_DIT_Encoder does for a complete input block.  A
    // one-coordinate shortened tail changes only the final four-row boundary;
    // the weighted primitive supplies mathematical zero without reading a
    // nonexistent public pointer and leaves the remaining traversal regular.
    for (uint32_t r = 0; r < side; r += 4)
    {
        if (!shortened || r + 4U <= data_count)
        {
            AVX2FF8IFFTButterfly4Out(
                data[r], data[r + 1], data[r + 2], data[r + 3],
                work[r], work[r + 1], work[r + 2], work[r + 3],
                inverse_skew[r + 1], inverse_skew[r + 3],
                inverse_skew[r + 2], byte_count);
        }
        else
        {
            LEO_DEBUG_ASSERT(r + 3U == data_count);
            AVX2FF8IFFTButterfly4OutLastZero(
                data[r], data[r + 1], data[r + 2],
                work[r], work[r + 1], work[r + 2], work[r + 3],
                inverse_skew[r + 1], inverse_skew[r + 3],
                inverse_skew[r + 2], byte_count);
        }
    }

    uint32_t distance = 4;
    uint32_t distance4 = 16;
    for (; distance4 <= side;
         distance = distance4, distance4 <<= 2)
    {
        for (uint32_t r = 0; r < side; r += distance4)
        {
            const uint32_t end = r + distance;
            AVX2FF8Butterfly4RangePrepared<true>(
                work + r, distance,
                inverse_skew[end], inverse_skew[end + distance * 2U],
                inverse_skew[end + distance], byte_count);
        }
    }
    if (distance < side)
    {
        AVX2FF8IFFTButterfly2RangePrepared(
            work, distance, inverse_skew[distance], byte_count);
    }

    // Evaluate the coefficient block on the legacy parity coset.  The
    // promoted production schedule uses the same fused radix-four operation
    // at every complete stage for these sizes.
    distance4 = side;
    distance = side >> 2;
    for (; distance != 0;
         distance4 = distance, distance >>= 2)
    {
        for (uint32_t r = 0; r < side; r += distance4)
        {
            const uint32_t end = r + distance;
            AVX2FF8Butterfly4RangePrepared<false>(
                work + r, distance,
                forward_skew[end], forward_skew[end + distance * 2U],
                forward_skew[end + distance], byte_count);
        }
    }
    if (distance4 == 2)
    {
        for (uint32_t r = 0; r < side; r += 2)
        {
            const uint16_t log = forward_skew[r + 1];
            if (log == kZeroSkew)
                AVX2XorMemory(work[r + 1], work[r], byte_count);
            else
                AVX2FF8Butterfly2<false>(
                    work[r], work[r + 1], log, byte_count);
        }
    }
}

#endif // LEO2_AVX512_VARIANT

#endif // LEO_HAS_FF8

#ifdef LEO_HAS_FF16
static LEO2_AVX2_FORCE_INLINE void AVX2FF16MultiplyAddPair(
    __m256i& destination_low,
    __m256i& destination_high,
    __m256i source_low,
    __m256i source_high,
    const FF16NibbleTable& table)
{
    const __m256i low_tables[4] = {
        BroadcastTable(table.low[0]), BroadcastTable(table.low[1]),
        BroadcastTable(table.low[2]), BroadcastTable(table.low[3])
    };
    const __m256i high_tables[4] = {
        BroadcastTable(table.high[0]), BroadcastTable(table.high[1]),
        BroadcastTable(table.high[2]), BroadcastTable(table.high[3])
    };
    __m256i product_low;
    __m256i product_high;
    AVX2FF16ProductVectors(source_low, source_high,
        low_tables, high_tables, product_low, product_high);
    destination_low = _mm256_xor_si256(destination_low, product_low);
    destination_high = _mm256_xor_si256(destination_high, product_high);
}

template<bool Inverse>
static void AVX2FF16Butterfly4Split(
    void* value0, void* value1, void* value2, void* value3,
    uint16_t log01, uint16_t log23, uint16_t log02,
    uint64_t byte_count)
{
    static const uint16_t kZeroSkew = 65535;
    if (Inverse)
    {
        if (log01 == kZeroSkew)
            AVX2XorMemory(value1, value0, byte_count);
        else
            AVX2FF16Butterfly2<true>(value0, value1, log01, byte_count);
        if (log23 == kZeroSkew)
            AVX2XorMemory(value3, value2, byte_count);
        else
            AVX2FF16Butterfly2<true>(value2, value3, log23, byte_count);
        if (log02 == kZeroSkew)
        {
            AVX2XorMemory(value2, value0, byte_count);
            AVX2XorMemory(value3, value1, byte_count);
        }
        else
        {
            AVX2FF16Butterfly2<true>(value0, value2, log02, byte_count);
            AVX2FF16Butterfly2<true>(value1, value3, log02, byte_count);
        }
    }
    else
    {
        if (log02 == kZeroSkew)
        {
            AVX2XorMemory(value2, value0, byte_count);
            AVX2XorMemory(value3, value1, byte_count);
        }
        else
        {
            AVX2FF16Butterfly2<false>(value0, value2, log02, byte_count);
            AVX2FF16Butterfly2<false>(value1, value3, log02, byte_count);
        }
        if (log01 == kZeroSkew)
            AVX2XorMemory(value1, value0, byte_count);
        else
            AVX2FF16Butterfly2<false>(value0, value1, log01, byte_count);
        if (log23 == kZeroSkew)
            AVX2XorMemory(value3, value2, byte_count);
        else
            AVX2FF16Butterfly2<false>(value2, value3, log23, byte_count);
    }
}

template<bool Inverse>
static void AVX2FF16Butterfly4(
    void* value0, void* value1, void* value2, void* value3,
    uint16_t log01, uint16_t log23, uint16_t log02,
    uint64_t byte_count)
{
    static const uint16_t kZeroSkew = 65535;
    // Keep small working sets fused.  This also covers public tails such as
    // 66 bytes, which staging rounds to 128 bytes.  Larger shards use the
    // split schedule to bound register pressure and table residency.
    if (byte_count < 64 || byte_count > 128)
    {
        AVX2FF16Butterfly4Split<Inverse>(
            value0, value1, value2, value3,
            log01, log23, log02, byte_count);
        return;
    }
    uint8_t* values[4] = {
        static_cast<uint8_t*>(value0), static_cast<uint8_t*>(value1),
        static_cast<uint8_t*>(value2), static_cast<uint8_t*>(value3)
    };
    uint64_t offset = 0;
    while (byte_count - offset >= 64)
    {
        __m256i low[4];
        __m256i high[4];
        for (unsigned lane = 0; lane < 4; ++lane)
        {
            low[lane] = _mm256_loadu_si256(
                reinterpret_cast<const __m256i*>(values[lane] + offset));
            high[lane] = _mm256_loadu_si256(
                reinterpret_cast<const __m256i*>(
                    values[lane] + offset + 32));
        }

        if (Inverse)
        {
            low[1] = _mm256_xor_si256(low[1], low[0]);
            high[1] = _mm256_xor_si256(high[1], high[0]);
            if (log01 != kZeroSkew)
                AVX2FF16MultiplyAddPair(low[0], high[0], low[1], high[1],
                    FF16Tables[log01]);

            low[3] = _mm256_xor_si256(low[3], low[2]);
            high[3] = _mm256_xor_si256(high[3], high[2]);
            if (log23 != kZeroSkew)
                AVX2FF16MultiplyAddPair(low[2], high[2], low[3], high[3],
                    FF16Tables[log23]);

            low[2] = _mm256_xor_si256(low[2], low[0]);
            high[2] = _mm256_xor_si256(high[2], high[0]);
            low[3] = _mm256_xor_si256(low[3], low[1]);
            high[3] = _mm256_xor_si256(high[3], high[1]);
            if (log02 != kZeroSkew)
            {
                AVX2FF16MultiplyAddPair(low[0], high[0], low[2], high[2],
                    FF16Tables[log02]);
                AVX2FF16MultiplyAddPair(low[1], high[1], low[3], high[3],
                    FF16Tables[log02]);
            }
        }
        else
        {
            if (log02 != kZeroSkew)
            {
                AVX2FF16MultiplyAddPair(low[0], high[0], low[2], high[2],
                    FF16Tables[log02]);
                AVX2FF16MultiplyAddPair(low[1], high[1], low[3], high[3],
                    FF16Tables[log02]);
            }
            low[2] = _mm256_xor_si256(low[2], low[0]);
            high[2] = _mm256_xor_si256(high[2], high[0]);
            low[3] = _mm256_xor_si256(low[3], low[1]);
            high[3] = _mm256_xor_si256(high[3], high[1]);

            if (log01 != kZeroSkew)
                AVX2FF16MultiplyAddPair(low[0], high[0], low[1], high[1],
                    FF16Tables[log01]);
            low[1] = _mm256_xor_si256(low[1], low[0]);
            high[1] = _mm256_xor_si256(high[1], high[0]);

            if (log23 != kZeroSkew)
                AVX2FF16MultiplyAddPair(low[2], high[2], low[3], high[3],
                    FF16Tables[log23]);
            low[3] = _mm256_xor_si256(low[3], low[2]);
            high[3] = _mm256_xor_si256(high[3], high[2]);
        }

        for (unsigned lane = 0; lane < 4; ++lane)
        {
            _mm256_storeu_si256(
                reinterpret_cast<__m256i*>(values[lane] + offset), low[lane]);
            _mm256_storeu_si256(reinterpret_cast<__m256i*>(
                values[lane] + offset + 32), high[lane]);
        }
        offset += 64;
    }

    const uint64_t residual = byte_count - offset;
    if (residual == 0)
        return;
    AVX2FF16Butterfly4Split<Inverse>(
        values[0] + offset, values[1] + offset,
        values[2] + offset, values[3] + offset,
        log01, log23, log02, residual);
}

static void AVX2FF16IFFTButterfly4(
    void* value0, void* value1, void* value2, void* value3,
    uint16_t log01, uint16_t log23, uint16_t log02,
    uint64_t byte_count)
{
    AVX2FF16Butterfly4<true>(value0, value1, value2, value3,
        log01, log23, log02, byte_count);
}

#undef LEO2_AVX2_FORCE_INLINE

static void AVX2FF16FFTButterfly4(
    void* value0, void* value1, void* value2, void* value3,
    uint16_t log01, uint16_t log23, uint16_t log02,
    uint64_t byte_count)
{
    AVX2FF16Butterfly4<false>(value0, value1, value2, value3,
        log01, log23, log02, byte_count);
}

template<bool Inverse>
static void AVX2FF16Butterfly4Out(
    const void* input0_pointer, const void* input1_pointer,
    const void* input2_pointer, const void* input3_pointer,
    void* output0_pointer, void* output1_pointer,
    void* output2_pointer, void* output3_pointer,
    uint16_t log01, uint16_t log23, uint16_t log02,
    uint64_t byte_count)
{
    static const uint16_t kZeroSkew = 65535;
    const uint8_t* inputs[4] = {
        static_cast<const uint8_t*>(input0_pointer),
        static_cast<const uint8_t*>(input1_pointer),
        static_cast<const uint8_t*>(input2_pointer),
        static_cast<const uint8_t*>(input3_pointer)
    };
    uint8_t* outputs[4] = {
        static_cast<uint8_t*>(output0_pointer),
        static_cast<uint8_t*>(output1_pointer),
        static_cast<uint8_t*>(output2_pointer),
        static_cast<uint8_t*>(output3_pointer)
    };
    uint64_t offset = 0;
    while (byte_count - offset >= 64)
    {
        __m256i low[4];
        __m256i high[4];
        for (unsigned lane = 0; lane < 4; ++lane)
        {
            low[lane] = _mm256_loadu_si256(
                reinterpret_cast<const __m256i*>(inputs[lane] + offset));
            high[lane] = _mm256_loadu_si256(
                reinterpret_cast<const __m256i*>(
                    inputs[lane] + offset + 32));
        }
        if (Inverse)
        {
            low[1] = _mm256_xor_si256(low[1], low[0]);
            high[1] = _mm256_xor_si256(high[1], high[0]);
            if (log01 != kZeroSkew)
                AVX2FF16MultiplyAddPair(low[0], high[0],
                    low[1], high[1], FF16Tables[log01]);
            low[3] = _mm256_xor_si256(low[3], low[2]);
            high[3] = _mm256_xor_si256(high[3], high[2]);
            if (log23 != kZeroSkew)
                AVX2FF16MultiplyAddPair(low[2], high[2],
                    low[3], high[3], FF16Tables[log23]);
            low[2] = _mm256_xor_si256(low[2], low[0]);
            high[2] = _mm256_xor_si256(high[2], high[0]);
            low[3] = _mm256_xor_si256(low[3], low[1]);
            high[3] = _mm256_xor_si256(high[3], high[1]);
            if (log02 != kZeroSkew)
            {
                AVX2FF16MultiplyAddPair(low[0], high[0],
                    low[2], high[2], FF16Tables[log02]);
                AVX2FF16MultiplyAddPair(low[1], high[1],
                    low[3], high[3], FF16Tables[log02]);
            }
        }
        else
        {
            if (log02 != kZeroSkew)
            {
                AVX2FF16MultiplyAddPair(low[0], high[0],
                    low[2], high[2], FF16Tables[log02]);
                AVX2FF16MultiplyAddPair(low[1], high[1],
                    low[3], high[3], FF16Tables[log02]);
            }
            low[2] = _mm256_xor_si256(low[2], low[0]);
            high[2] = _mm256_xor_si256(high[2], high[0]);
            low[3] = _mm256_xor_si256(low[3], low[1]);
            high[3] = _mm256_xor_si256(high[3], high[1]);
            if (log01 != kZeroSkew)
                AVX2FF16MultiplyAddPair(low[0], high[0],
                    low[1], high[1], FF16Tables[log01]);
            low[1] = _mm256_xor_si256(low[1], low[0]);
            high[1] = _mm256_xor_si256(high[1], high[0]);
            if (log23 != kZeroSkew)
                AVX2FF16MultiplyAddPair(low[2], high[2],
                    low[3], high[3], FF16Tables[log23]);
            low[3] = _mm256_xor_si256(low[3], low[2]);
            high[3] = _mm256_xor_si256(high[3], high[2]);
        }
        for (unsigned lane = 0; lane < 4; ++lane)
        {
            _mm256_storeu_si256(
                reinterpret_cast<__m256i*>(outputs[lane] + offset), low[lane]);
            _mm256_storeu_si256(reinterpret_cast<__m256i*>(
                outputs[lane] + offset + 32), high[lane]);
        }
        offset += 64;
    }
    const uint64_t symbols = (byte_count - offset) / 2;
    for (uint64_t i = 0; i < symbols; ++i)
    {
        uint16_t x[4];
        for (unsigned lane = 0; lane < 4; ++lane)
            x[lane] = static_cast<uint16_t>(inputs[lane][offset + i] |
                (static_cast<unsigned>(inputs[lane][offset + symbols + i]) << 8));
        if (Inverse)
        {
            x[1] ^= x[0];
            if (log01 != kZeroSkew)
                x[0] ^= FF16Product(log01, x[1]);
            x[3] ^= x[2];
            if (log23 != kZeroSkew)
                x[2] ^= FF16Product(log23, x[3]);
            x[2] ^= x[0];
            x[3] ^= x[1];
            if (log02 != kZeroSkew)
            {
                x[0] ^= FF16Product(log02, x[2]);
                x[1] ^= FF16Product(log02, x[3]);
            }
        }
        else
        {
            if (log02 != kZeroSkew)
            {
                x[0] ^= FF16Product(log02, x[2]);
                x[1] ^= FF16Product(log02, x[3]);
            }
            x[2] ^= x[0];
            x[3] ^= x[1];
            if (log01 != kZeroSkew)
                x[0] ^= FF16Product(log01, x[1]);
            x[1] ^= x[0];
            if (log23 != kZeroSkew)
                x[2] ^= FF16Product(log23, x[3]);
            x[3] ^= x[2];
        }
        for (unsigned lane = 0; lane < 4; ++lane)
        {
            outputs[lane][offset + i] = static_cast<uint8_t>(x[lane]);
            outputs[lane][offset + symbols + i] =
                static_cast<uint8_t>(x[lane] >> 8);
        }
    }
}

static void AVX2FF16IFFTButterfly4Out(
    const void* input0, const void* input1,
    const void* input2, const void* input3,
    void* output0, void* output1, void* output2, void* output3,
    uint16_t log01, uint16_t log23, uint16_t log02,
    uint64_t byte_count)
{
    AVX2FF16Butterfly4Out<true>(
        input0, input1, input2, input3,
        output0, output1, output2, output3,
        log01, log23, log02, byte_count);
}

static void AVX2FF16FFTButterfly4Out(
    const void* input0, const void* input1,
    const void* input2, const void* input3,
    void* output0, void* output1, void* output2, void* output3,
    uint16_t log01, uint16_t log23, uint16_t log02,
    uint64_t byte_count)
{
    AVX2FF16Butterfly4Out<false>(
        input0, input1, input2, input3,
        output0, output1, output2, output3,
        log01, log23, log02, byte_count);
}

#if !defined(LEO2_AVX512_VARIANT)
template<bool Inverse>
static void AVX2FF16Butterfly2RangePrepared(
    void* const* work,
    unsigned x_base,
    unsigned y_base,
    unsigned distance,
    uint16_t multiplier_log,
    uint64_t byte_count)
{
    // Every pair in a transform range uses the same fixed multiplier.  Keep
    // its eight nibble rows prepared across the independent pairs instead of
    // broadcasting them again for every shard pair.  The AVX2 product still
    // uses the established byte order and arithmetic within each butterfly.
    static const uint16_t kZeroSkew = 65535;
    if (multiplier_log == kZeroSkew)
    {
        for (unsigned i = 0; i < distance; ++i)
            AVX2XorMemory(work[y_base + i], work[x_base + i], byte_count);
        return;
    }
    const FF16NibbleTable& table = FF16Tables[multiplier_log];
    const __m256i low_tables[4] = {
        BroadcastTable(table.low[0]), BroadcastTable(table.low[1]),
        BroadcastTable(table.low[2]), BroadcastTable(table.low[3])
    };
    const __m256i high_tables[4] = {
        BroadcastTable(table.high[0]), BroadcastTable(table.high[1]),
        BroadcastTable(table.high[2]), BroadcastTable(table.high[3])
    };
    for (unsigned i = 0; i < distance; ++i)
    {
        AVX2FF16Butterfly2Prepared<Inverse>(
            work[x_base + i], work[y_base + i], multiplier_log,
            low_tables, high_tables, byte_count);
    }
}
#endif

template<bool Inverse>
static void AVX2FF16Butterfly4Range(
    void* const* work, unsigned distance,
    uint16_t log01, uint16_t log23, uint16_t log02,
    uint64_t byte_count, bool prefer_fused)
{
#if defined(LEO2_AVX512_VARIANT)
    for (unsigned i = 0; i < distance; ++i)
    {
        void* const value0 = work[i];
        void* const value1 = work[i + distance];
        void* const value2 = work[i + distance * 2U];
        void* const value3 = work[i + distance * 3U];
        if (prefer_fused)
            AVX2FF16Butterfly4<Inverse>(value0, value1, value2, value3,
                log01, log23, log02, byte_count);
        else
            AVX2FF16Butterfly4Split<Inverse>(
                value0, value1, value2, value3,
                log01, log23, log02, byte_count);
    }
#else
    // The Butterfly4Range contract makes all coordinates disjoint, so the
    // split radix-four layers are independent across i.  Execute each layer
    // as a range so its fixed tables remain prepared, then preserve the
    // inverse/forward layer order for every coordinate.  AVX-512VL retains
    // the fused/per-pair schedule above: its expanded register file has a
    // different measured crossover and does not need this AVX2 workaround.
    if (prefer_fused)
    {
        for (unsigned i = 0; i < distance; ++i)
        {
            AVX2FF16Butterfly4<Inverse>(
                work[i], work[i + distance],
                work[i + distance * 2U], work[i + distance * 3U],
                log01, log23, log02, byte_count);
        }
        return;
    }
    if (Inverse)
    {
        AVX2FF16Butterfly2RangePrepared<true>(
            work, 0, distance, distance, log01, byte_count);
        AVX2FF16Butterfly2RangePrepared<true>(
            work, distance * 2U, distance * 3U,
            distance, log23, byte_count);
        AVX2FF16Butterfly2RangePrepared<true>(
            work, 0, distance * 2U, distance, log02, byte_count);
        AVX2FF16Butterfly2RangePrepared<true>(
            work, distance, distance * 3U,
            distance, log02, byte_count);
    }
    else
    {
        AVX2FF16Butterfly2RangePrepared<false>(
            work, 0, distance * 2U, distance, log02, byte_count);
        AVX2FF16Butterfly2RangePrepared<false>(
            work, distance, distance * 3U,
            distance, log02, byte_count);
        AVX2FF16Butterfly2RangePrepared<false>(
            work, 0, distance, distance, log01, byte_count);
        AVX2FF16Butterfly2RangePrepared<false>(
            work, distance * 2U, distance * 3U,
            distance, log23, byte_count);
    }
#endif
}

static void AVX2FF16IFFTButterfly4Range(
    void* const* work, unsigned distance,
    uint16_t log01, uint16_t log23, uint16_t log02,
    uint64_t byte_count, bool prefer_fused)
{
    AVX2FF16Butterfly4Range<true>(work, distance,
        log01, log23, log02, byte_count, prefer_fused);
}

static void AVX2FF16FFTButterfly4Range(
    void* const* work, unsigned distance,
    uint16_t log01, uint16_t log23, uint16_t log02,
    uint64_t byte_count, bool prefer_fused)
{
    AVX2FF16Butterfly4Range<false>(work, distance,
        log01, log23, log02, byte_count, prefer_fused);
}
#endif // LEO_HAS_FF16

static const Ops AVX2Ops = {
    LEO2_AVX_BACKEND_KIND,
    LEO2_AVX_BACKEND_NAME,
#ifdef LEO_HAS_FF8
    AVX2FF8Multiply,
    AVX2FF8MultiplyAdd,
#else
    NULL,
    NULL,
#endif
#ifdef LEO_HAS_FF16
    AVX2FF16Multiply,
    AVX2FF16MultiplyAdd,
#else
    NULL,
    NULL,
#endif
    AVX2XorMemory,
    AVX2XorMemory2To1,
    AVX2XorMemorySources,
    AVX2XorMemory4,
#ifdef LEO_HAS_FF8
    AVX2FF8IFFTButterfly2,
    AVX2FF8FFTButterfly2,
    AVX2FF8FFTButterfly2Out,
    AVX2FF8IFFTButterfly2Xor,
    AVX2FF8IFFTButterfly4,
    AVX2FF8FFTButterfly4,
    AVX2FF8IFFTButterfly4Out,
    AVX2FF8WeightedIFFTButterfly4,
    AVX2FF8FFTButterfly4Out,
    AVX2FF8IFFTButterfly4Range,
    AVX2FF8FFTButterfly4Range,
    AVX2FF8IFFTButterfly4XorRange,
#else
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
#endif
#ifdef LEO_HAS_FF16
    AVX2FF16IFFTButterfly2,
    AVX2FF16FFTButterfly2,
    AVX2FF16FFTButterfly2Out,
    AVX2FF16IFFTButterfly2Xor,
    AVX2FF16IFFTButterfly4,
    AVX2FF16FFTButterfly4,
    AVX2FF16IFFTButterfly4Out,
    AVX2FF16FFTButterfly4Out,
    AVX2FF16IFFTButterfly4Range,
    AVX2FF16FFTButterfly4Range
#else
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL
#endif
#if defined(LEO_HAS_FF8) && defined(LEO2_AVX512_VARIANT)
    , AVX2FF8HighEncodeOneBlock
#else
    , NULL
#endif
#if defined(LEO_HAS_FF8) && !defined(LEO2_AVX512_VARIANT)
    , AVX2FF8HighEncodeSmall
#else
    , NULL
#endif
};

const Ops* LEO2_AVX_INITIALIZER(const InitializeArgs& args)
{
#ifdef LEO_HAS_FF8
    if (!args.ff8_multiply_log)
        return NULL;
#endif
#ifdef LEO_HAS_FF16
    if (!args.ff16_multiply_log)
        return NULL;
#endif
#if defined(LEO2_AVX512_VARIANT)
# ifdef LEO2_ENABLE_TEST_HOOKS
#  ifdef LEO_HAS_FF8
    if (TestShouldFailAllocation(LEO2_AVX_BACKEND_KIND, false))
        return NULL;
#  endif
#  ifdef LEO_HAS_FF16
    if (TestShouldFailAllocation(LEO2_AVX_BACKEND_KIND, true))
        return NULL;
#  endif
# endif
    const Ops* const base = InitializeAVX2(args);
    if (!base)
        return NULL;
# ifdef LEO_HAS_FF8
    FF8Tables = static_cast<const FF8NibbleTable*>(GetAVX2FF8Tables());
    if (!FF8Tables)
        return NULL;
# endif
# ifdef LEO_HAS_FF16
    FF16Tables = static_cast<const FF16NibbleTable*>(GetAVX2FF16Tables());
    if (!FF16Tables)
        return NULL;
# endif
    return &AVX2Ops;
#else
#if defined(LEO_HAS_FF8) && defined(LEO_HAS_FF16)
    if (FF8Tables || FF16Tables)
        return FF8Tables && FF16Tables ? &AVX2Ops : NULL;
#elif defined(LEO_HAS_FF8)
    if (FF8Tables)
        return &AVX2Ops;
#else
    if (FF16Tables)
        return &AVX2Ops;
#endif
#ifdef LEO_HAS_FF8
#ifdef LEO2_ENABLE_TEST_HOOKS
    if (TestShouldFailAllocation(LEO2_AVX_BACKEND_KIND, false))
        return NULL;
#endif
    std::unique_ptr<FF8NibbleTable[]> ff8(
        new (std::nothrow) FF8NibbleTable[256]);
    if (!ff8)
        return NULL;
#endif
#ifdef LEO_HAS_FF16
#ifdef LEO2_ENABLE_TEST_HOOKS
    if (TestShouldFailAllocation(LEO2_AVX_BACKEND_KIND, true))
        return NULL;
#endif
    std::unique_ptr<FF16NibbleTable[]> ff16(
        new (std::nothrow) FF16NibbleTable[65536]);
    if (!ff16)
        return NULL;
#endif

#ifdef LEO_HAS_FF8
    for (unsigned log = 0; log < 256; ++log)
    {
        for (unsigned value = 0; value < 16; ++value)
        {
            ff8[log].low[value] = args.ff8_multiply_log(
                static_cast<uint8_t>(value), static_cast<uint8_t>(log));
            ff8[log].high[value] = args.ff8_multiply_log(
                static_cast<uint8_t>(value << 4), static_cast<uint8_t>(log));
        }
    }
#endif

#ifdef LEO_HAS_FF16
    for (int log = 0; log < 65536; ++log)
    {
        for (unsigned nibble = 0; nibble < 4; ++nibble)
        {
            for (unsigned value = 0; value < 16; ++value)
            {
                const uint16_t product = args.ff16_multiply_log(
                    static_cast<uint16_t>(value << (nibble * 4)),
                    static_cast<uint16_t>(log));
                ff16[log].low[nibble][value] =
                    static_cast<uint8_t>(product);
                ff16[log].high[nibble][value] =
                    static_cast<uint8_t>(product >> 8);
            }
        }
    }
#endif
#ifdef LEO_HAS_FF8
    FF8Tables = ff8.release();
#endif
#ifdef LEO_HAS_FF16
    FF16Tables = ff16.release();
#endif
    return &AVX2Ops;
#endif
}

#if !defined(LEO2_AVX512_VARIANT)
const void* GetAVX2FF8Tables()
{
#ifdef LEO_HAS_FF8
    return FF8Tables;
#else
    return NULL;
#endif
}

const void* GetAVX2FF16Tables()
{
#ifdef LEO_HAS_FF16
    return FF16Tables;
#else
    return NULL;
#endif
}
#endif

#ifdef LEO2_ENABLE_TEST_HOOKS
void LEO2_AVX_TABLE_STATE(TestBackendState* state)
{
#ifdef LEO_HAS_FF8
    state->ff8_published = FF8Tables != NULL;
    state->ff8_bytes = 256U * sizeof(FF8NibbleTable);
#else
    state->ff8_published = false;
    state->ff8_bytes = 0;
#endif
#ifdef LEO_HAS_FF16
    state->ff16_published = FF16Tables != NULL;
    state->ff16_bytes = 65536U * sizeof(FF16NibbleTable);
#else
    state->ff16_published = false;
    state->ff16_bytes = 0;
#endif
}
#endif

#undef LEO2_AVX_TABLE_STATE
#undef LEO2_AVX_INITIALIZER
#undef LEO2_AVX_BACKEND_NAME
#undef LEO2_AVX_BACKEND_KIND

}} // namespace leopard::backend
