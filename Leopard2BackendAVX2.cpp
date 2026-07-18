/*
    Copyright (c) 2017 Christopher A. Taylor.  All rights reserved.
    See Leopard2Backend.h for the full BSD license notice.
*/

#include "Leopard2Backend.h"
#include "LeopardCommon.h"

#include <immintrin.h>
#include <memory>
#include <new>

namespace leopard { namespace backend {

#ifdef LEO_HAS_FF8
struct FF8NibbleTable
{
    uint8_t low[16];
    uint8_t high[16];
};
static FF8NibbleTable* FF8Tables = NULL;
#endif

#ifdef LEO_HAS_FF16
struct FF16NibbleTable
{
    uint8_t low[4][16];
    uint8_t high[4][16];
};
static FF16NibbleTable* FF16Tables = NULL;
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
static void AVX2FF16Butterfly2(
    void* x_pointer,
    void* y_pointer,
    uint16_t multiplier_log,
    uint64_t byte_count)
{
    uint8_t* x = static_cast<uint8_t*>(x_pointer);
    uint8_t* y = static_cast<uint8_t*>(y_pointer);
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

static void AVX2FF8FFTButterfly4Out(
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
        __m256i x3 = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(input3));
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
        _mm256_storeu_si256(reinterpret_cast<__m256i*>(output0), x0);
        _mm256_storeu_si256(reinterpret_cast<__m256i*>(output1), x1);
        _mm256_storeu_si256(reinterpret_cast<__m256i*>(output2), x2);
        _mm256_storeu_si256(reinterpret_cast<__m256i*>(output3), x3);
        input0 += 32;
        input1 += 32;
        input2 += 32;
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
        uint8_t x3 = *input3++;
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
        *output0++ = x0;
        *output1++ = x1;
        *output2++ = x2;
        *output3++ = x3;
    }
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
    (void)prefer_fused;
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

static void AVX2FF16FFTButterfly4Out(
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
        for (unsigned lane = 0; lane < 4; ++lane)
        {
            outputs[lane][offset + i] = static_cast<uint8_t>(x[lane]);
            outputs[lane][offset + symbols + i] =
                static_cast<uint8_t>(x[lane] >> 8);
        }
    }
}

template<bool Inverse>
static void AVX2FF16Butterfly4Range(
    void* const* work, unsigned distance,
    uint16_t log01, uint16_t log23, uint16_t log02,
    uint64_t byte_count, bool prefer_fused)
{
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
    LEO2_BACKEND_AVX2,
    "avx2",
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
    AVX2XorMemory4,
#ifdef LEO_HAS_FF8
    AVX2FF8IFFTButterfly2,
    AVX2FF8FFTButterfly2,
    AVX2FF8FFTButterfly2Out,
    AVX2FF8IFFTButterfly2Xor,
    AVX2FF8IFFTButterfly4,
    AVX2FF8FFTButterfly4,
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
#endif
#ifdef LEO_HAS_FF16
    AVX2FF16IFFTButterfly2,
    AVX2FF16FFTButterfly2,
    AVX2FF16FFTButterfly2Out,
    AVX2FF16IFFTButterfly4,
    AVX2FF16FFTButterfly4,
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
    NULL
#endif
};

const Ops* InitializeAVX2(const InitializeArgs& args)
{
#ifdef LEO_HAS_FF8
    if (!args.ff8_multiply_log)
        return NULL;
#endif
#ifdef LEO_HAS_FF16
    if (!args.ff16_multiply_log)
        return NULL;
#endif
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
    if (TestShouldFailAllocation(LEO2_BACKEND_AVX2, false))
        return NULL;
#endif
    std::unique_ptr<FF8NibbleTable[]> ff8(
        new (std::nothrow) FF8NibbleTable[256]);
    if (!ff8)
        return NULL;
#endif
#ifdef LEO_HAS_FF16
#ifdef LEO2_ENABLE_TEST_HOOKS
    if (TestShouldFailAllocation(LEO2_BACKEND_AVX2, true))
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
}

#ifdef LEO2_ENABLE_TEST_HOOKS
void TestGetAVX2TableState(TestBackendState* state)
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

}} // namespace leopard::backend
