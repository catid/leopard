/*
    Copyright (c) 2017 Christopher A. Taylor.  All rights reserved.
    See Leopard2Backend.h for the full BSD license notice.
*/

#include "Leopard2Backend.h"

#include <new>
#include <tmmintrin.h>

namespace leopard { namespace backend {

struct FF8NibbleTable
{
    uint8_t low[16];
    uint8_t high[16];
};

struct FF16NibbleTable
{
    uint8_t low[4][16];
    uint8_t high[4][16];
};

static FF8NibbleTable* FF8Tables = NULL;
static FF16NibbleTable* FF16Tables = NULL;

static uint8_t FF8Product(uint16_t log, uint8_t value)
{
    const FF8NibbleTable& table = FF8Tables[log];
    return static_cast<uint8_t>(
        table.low[value & 15U] ^ table.high[value >> 4]);
}

template<bool Add>
static void SSSE3FF8Operation(
    void* destination,
    const void* source,
    uint16_t multiplier_log,
    uint64_t byte_count)
{
    uint8_t* output = static_cast<uint8_t*>(destination);
    const uint8_t* input = static_cast<const uint8_t*>(source);
    const FF8NibbleTable& table = FF8Tables[multiplier_log];
    const __m128i low_table = _mm_loadu_si128(
        reinterpret_cast<const __m128i*>(table.low));
    const __m128i high_table = _mm_loadu_si128(
        reinterpret_cast<const __m128i*>(table.high));
    const __m128i nibble_mask = _mm_set1_epi8(15);
    while (byte_count >= 16)
    {
        const __m128i data = _mm_loadu_si128(
            reinterpret_cast<const __m128i*>(input));
        const __m128i low = _mm_shuffle_epi8(
            low_table, _mm_and_si128(data, nibble_mask));
        const __m128i high = _mm_shuffle_epi8(high_table,
            _mm_and_si128(_mm_srli_epi64(data, 4), nibble_mask));
        __m128i product = _mm_xor_si128(low, high);
        if (Add)
            product = _mm_xor_si128(product, _mm_loadu_si128(
                reinterpret_cast<const __m128i*>(output)));
        _mm_storeu_si128(reinterpret_cast<__m128i*>(output), product);
        input += 16;
        output += 16;
        byte_count -= 16;
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

static void SSSE3FF8Multiply(
    void* destination, const void* source,
    uint16_t multiplier_log, uint64_t byte_count)
{
    SSSE3FF8Operation<false>(
        destination, source, multiplier_log, byte_count);
}

static void SSSE3FF8MultiplyAdd(
    void* destination, const void* source,
    uint16_t multiplier_log, uint64_t byte_count)
{
    SSSE3FF8Operation<true>(
        destination, source, multiplier_log, byte_count);
}

static __m128i SSSE3FF8ProductVector(
    __m128i data,
    __m128i low_table,
    __m128i high_table)
{
    const __m128i nibble_mask = _mm_set1_epi8(15);
    const __m128i low = _mm_shuffle_epi8(
        low_table, _mm_and_si128(data, nibble_mask));
    const __m128i high = _mm_shuffle_epi8(high_table,
        _mm_and_si128(_mm_srli_epi64(data, 4), nibble_mask));
    return _mm_xor_si128(low, high);
}

template<bool Inverse>
static void SSSE3FF8Butterfly2(
    void* x_pointer,
    void* y_pointer,
    uint16_t multiplier_log,
    uint64_t byte_count)
{
    uint8_t* x = static_cast<uint8_t*>(x_pointer);
    uint8_t* y = static_cast<uint8_t*>(y_pointer);
    const FF8NibbleTable& table = FF8Tables[multiplier_log];
    const __m128i low_table = _mm_loadu_si128(
        reinterpret_cast<const __m128i*>(table.low));
    const __m128i high_table = _mm_loadu_si128(
        reinterpret_cast<const __m128i*>(table.high));
    while (byte_count >= 16)
    {
        __m128i x_value = _mm_loadu_si128(
            reinterpret_cast<const __m128i*>(x));
        __m128i y_value = _mm_loadu_si128(
            reinterpret_cast<const __m128i*>(y));
        if (Inverse)
        {
            y_value = _mm_xor_si128(y_value, x_value);
            x_value = _mm_xor_si128(x_value,
                SSSE3FF8ProductVector(y_value, low_table, high_table));
        }
        else
        {
            x_value = _mm_xor_si128(x_value,
                SSSE3FF8ProductVector(y_value, low_table, high_table));
            y_value = _mm_xor_si128(y_value, x_value);
        }
        _mm_storeu_si128(reinterpret_cast<__m128i*>(x), x_value);
        _mm_storeu_si128(reinterpret_cast<__m128i*>(y), y_value);
        x += 16;
        y += 16;
        byte_count -= 16;
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

static void SSSE3FF8IFFTButterfly2(
    void* x, void* y, uint16_t multiplier_log, uint64_t byte_count)
{
    SSSE3FF8Butterfly2<true>(x, y, multiplier_log, byte_count);
}

static void SSSE3FF8FFTButterfly2(
    void* x, void* y, uint16_t multiplier_log, uint64_t byte_count)
{
    SSSE3FF8Butterfly2<false>(x, y, multiplier_log, byte_count);
}

static void SSSE3FF8IFFTButterfly2Xor(
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
    const __m128i low_table = _mm_loadu_si128(
        reinterpret_cast<const __m128i*>(table.low));
    const __m128i high_table = _mm_loadu_si128(
        reinterpret_cast<const __m128i*>(table.high));
    while (byte_count >= 16)
    {
        const __m128i x_original = _mm_loadu_si128(
            reinterpret_cast<const __m128i*>(x_input));
        const __m128i y_value = _mm_xor_si128(x_original,
            _mm_loadu_si128(reinterpret_cast<const __m128i*>(y_input)));
        const __m128i x_value = _mm_xor_si128(x_original,
            SSSE3FF8ProductVector(y_value, low_table, high_table));
        const __m128i x_result = _mm_xor_si128(x_value,
            _mm_loadu_si128(reinterpret_cast<const __m128i*>(x_output)));
        const __m128i y_result = _mm_xor_si128(y_value,
            _mm_loadu_si128(reinterpret_cast<const __m128i*>(y_output)));
        _mm_storeu_si128(reinterpret_cast<__m128i*>(x_output), x_result);
        _mm_storeu_si128(reinterpret_cast<__m128i*>(y_output), y_result);
        x_input += 16;
        y_input += 16;
        x_output += 16;
        y_output += 16;
        byte_count -= 16;
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

static void SSSE3FF16ProductVectors(
    __m128i low_data,
    __m128i high_data,
    const FF16NibbleTable& table,
    __m128i& product_low,
    __m128i& product_high)
{
    const __m128i mask = _mm_set1_epi8(15);
    const __m128i nibbles[4] = {
        _mm_and_si128(low_data, mask),
        _mm_and_si128(_mm_srli_epi64(low_data, 4), mask),
        _mm_and_si128(high_data, mask),
        _mm_and_si128(_mm_srli_epi64(high_data, 4), mask)
    };
    product_low = _mm_setzero_si128();
    product_high = _mm_setzero_si128();
    for (unsigned i = 0; i < 4; ++i)
    {
        const __m128i low_table = _mm_loadu_si128(
            reinterpret_cast<const __m128i*>(table.low[i]));
        const __m128i high_table = _mm_loadu_si128(
            reinterpret_cast<const __m128i*>(table.high[i]));
        product_low = _mm_xor_si128(product_low,
            _mm_shuffle_epi8(low_table, nibbles[i]));
        product_high = _mm_xor_si128(product_high,
            _mm_shuffle_epi8(high_table, nibbles[i]));
    }
}

static void SSSE3FF16Block(
    uint8_t* output_low,
    uint8_t* output_high,
    const uint8_t* input_low,
    const uint8_t* input_high,
    const FF16NibbleTable& table,
    bool add)
{
    const __m128i mask = _mm_set1_epi8(15);
    const __m128i low_data = _mm_loadu_si128(
        reinterpret_cast<const __m128i*>(input_low));
    const __m128i high_data = _mm_loadu_si128(
        reinterpret_cast<const __m128i*>(input_high));
    const __m128i nibbles[4] = {
        _mm_and_si128(low_data, mask),
        _mm_and_si128(_mm_srli_epi64(low_data, 4), mask),
        _mm_and_si128(high_data, mask),
        _mm_and_si128(_mm_srli_epi64(high_data, 4), mask)
    };
    __m128i product_low = _mm_setzero_si128();
    __m128i product_high = _mm_setzero_si128();
    for (unsigned i = 0; i < 4; ++i)
    {
        const __m128i low_table = _mm_loadu_si128(
            reinterpret_cast<const __m128i*>(table.low[i]));
        const __m128i high_table = _mm_loadu_si128(
            reinterpret_cast<const __m128i*>(table.high[i]));
        product_low = _mm_xor_si128(product_low,
            _mm_shuffle_epi8(low_table, nibbles[i]));
        product_high = _mm_xor_si128(product_high,
            _mm_shuffle_epi8(high_table, nibbles[i]));
    }
    if (add)
    {
        product_low = _mm_xor_si128(product_low, _mm_loadu_si128(
            reinterpret_cast<const __m128i*>(output_low)));
        product_high = _mm_xor_si128(product_high, _mm_loadu_si128(
            reinterpret_cast<const __m128i*>(output_high)));
    }
    _mm_storeu_si128(reinterpret_cast<__m128i*>(output_low), product_low);
    _mm_storeu_si128(reinterpret_cast<__m128i*>(output_high), product_high);
}

template<bool Add>
static void SSSE3FF16Operation(
    void* destination,
    const void* source,
    uint16_t multiplier_log,
    uint64_t byte_count)
{
    uint8_t* output = static_cast<uint8_t*>(destination);
    const uint8_t* input = static_cast<const uint8_t*>(source);
    const FF16NibbleTable& table = FF16Tables[multiplier_log];
    uint64_t offset = 0;
    while (byte_count - offset >= 64)
    {
        SSSE3FF16Block(
            output + offset, output + offset + 32,
            input + offset, input + offset + 32, table, Add);
        SSSE3FF16Block(
            output + offset + 16, output + offset + 48,
            input + offset + 16, input + offset + 48, table, Add);
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

static void SSSE3FF16Multiply(
    void* destination, const void* source,
    uint16_t multiplier_log, uint64_t byte_count)
{
    SSSE3FF16Operation<false>(
        destination, source, multiplier_log, byte_count);
}

static void SSSE3FF16MultiplyAdd(
    void* destination, const void* source,
    uint16_t multiplier_log, uint64_t byte_count)
{
    SSSE3FF16Operation<true>(
        destination, source, multiplier_log, byte_count);
}

template<bool Inverse>
static void SSSE3FF16ButterflyBlock(
    uint8_t* x_low,
    uint8_t* x_high,
    uint8_t* y_low,
    uint8_t* y_high,
    const FF16NibbleTable& table)
{
    __m128i x_low_value = _mm_loadu_si128(
        reinterpret_cast<const __m128i*>(x_low));
    __m128i x_high_value = _mm_loadu_si128(
        reinterpret_cast<const __m128i*>(x_high));
    __m128i y_low_value = _mm_loadu_si128(
        reinterpret_cast<const __m128i*>(y_low));
    __m128i y_high_value = _mm_loadu_si128(
        reinterpret_cast<const __m128i*>(y_high));
    if (Inverse)
    {
        y_low_value = _mm_xor_si128(y_low_value, x_low_value);
        y_high_value = _mm_xor_si128(y_high_value, x_high_value);
    }
    __m128i product_low;
    __m128i product_high;
    SSSE3FF16ProductVectors(
        y_low_value, y_high_value, table, product_low, product_high);
    x_low_value = _mm_xor_si128(x_low_value, product_low);
    x_high_value = _mm_xor_si128(x_high_value, product_high);
    if (!Inverse)
    {
        y_low_value = _mm_xor_si128(y_low_value, x_low_value);
        y_high_value = _mm_xor_si128(y_high_value, x_high_value);
    }
    _mm_storeu_si128(reinterpret_cast<__m128i*>(x_low), x_low_value);
    _mm_storeu_si128(reinterpret_cast<__m128i*>(x_high), x_high_value);
    _mm_storeu_si128(reinterpret_cast<__m128i*>(y_low), y_low_value);
    _mm_storeu_si128(reinterpret_cast<__m128i*>(y_high), y_high_value);
}

template<bool Inverse>
static void SSSE3FF16Butterfly2(
    void* x_pointer,
    void* y_pointer,
    uint16_t multiplier_log,
    uint64_t byte_count)
{
    uint8_t* x = static_cast<uint8_t*>(x_pointer);
    uint8_t* y = static_cast<uint8_t*>(y_pointer);
    const FF16NibbleTable& table = FF16Tables[multiplier_log];
    uint64_t offset = 0;
    while (byte_count - offset >= 64)
    {
        SSSE3FF16ButterflyBlock<Inverse>(
            x + offset, x + offset + 32,
            y + offset, y + offset + 32, table);
        SSSE3FF16ButterflyBlock<Inverse>(
            x + offset + 16, x + offset + 48,
            y + offset + 16, y + offset + 48, table);
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

static void SSSE3FF16IFFTButterfly2(
    void* x, void* y, uint16_t multiplier_log, uint64_t byte_count)
{
    SSSE3FF16Butterfly2<true>(x, y, multiplier_log, byte_count);
}

static void SSSE3FF16FFTButterfly2(
    void* x, void* y, uint16_t multiplier_log, uint64_t byte_count)
{
    SSSE3FF16Butterfly2<false>(x, y, multiplier_log, byte_count);
}

static void SSSE3XorMemory(
    void* destination, const void* source, uint64_t byte_count)
{
    uint8_t* output = static_cast<uint8_t*>(destination);
    const uint8_t* input = static_cast<const uint8_t*>(source);
    while (byte_count >= 16)
    {
        const __m128i result = _mm_xor_si128(
            _mm_loadu_si128(reinterpret_cast<const __m128i*>(output)),
            _mm_loadu_si128(reinterpret_cast<const __m128i*>(input)));
        _mm_storeu_si128(reinterpret_cast<__m128i*>(output), result);
        output += 16;
        input += 16;
        byte_count -= 16;
    }
    while (byte_count-- != 0)
        *output++ ^= *input++;
}

static const Ops SSSE3Ops = {
    LEO2_BACKEND_SSSE3,
    "ssse3",
    SSSE3FF8Multiply,
    SSSE3FF8MultiplyAdd,
    SSSE3FF16Multiply,
    SSSE3FF16MultiplyAdd,
    SSSE3XorMemory,
    SSSE3FF8IFFTButterfly2,
    SSSE3FF8FFTButterfly2,
    SSSE3FF8IFFTButterfly2Xor,
    SSSE3FF16IFFTButterfly2,
    SSSE3FF16FFTButterfly2
};

const Ops* InitializeSSSE3(const InitializeArgs& args)
{
    if (!args.ff8_multiply_log || !args.ff16_multiply_log)
        return NULL;
    if (!FF8Tables)
        FF8Tables = new (std::nothrow) FF8NibbleTable[256];
    if (!FF16Tables)
        FF16Tables = new (std::nothrow) FF16NibbleTable[65536];
    if (!FF8Tables || !FF16Tables)
        return NULL;

    for (unsigned log = 0; log < 256; ++log)
    {
        for (unsigned value = 0; value < 16; ++value)
        {
            FF8Tables[log].low[value] = args.ff8_multiply_log(
                static_cast<uint8_t>(value), static_cast<uint8_t>(log));
            FF8Tables[log].high[value] = args.ff8_multiply_log(
                static_cast<uint8_t>(value << 4), static_cast<uint8_t>(log));
        }
    }

    for (int log = 0; log < 65536; ++log)
    {
        for (unsigned nibble = 0; nibble < 4; ++nibble)
        {
            for (unsigned value = 0; value < 16; ++value)
            {
                const uint16_t product = args.ff16_multiply_log(
                    static_cast<uint16_t>(value << (nibble * 4)),
                    static_cast<uint16_t>(log));
                FF16Tables[log].low[nibble][value] =
                    static_cast<uint8_t>(product);
                FF16Tables[log].high[nibble][value] =
                    static_cast<uint8_t>(product >> 8);
            }
        }
    }
    return &SSSE3Ops;
}

}} // namespace leopard::backend
