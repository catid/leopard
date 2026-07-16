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
    SSSE3XorMemory
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
