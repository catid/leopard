/*
    Copyright (c) 2017 Christopher A. Taylor.  All rights reserved.
    See Leopard2Backend.h for the full BSD license notice.
*/

#include "Leopard2Backend.h"

#include <immintrin.h>
#include <new>

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

static __m256i BroadcastTable(const uint8_t table[16])
{
    return _mm256_broadcastsi128_si256(_mm_loadu_si128(
        reinterpret_cast<const __m128i*>(table)));
}

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

static const Ops AVX2Ops = {
    LEO2_BACKEND_AVX2,
    "avx2",
    AVX2FF8Multiply,
    AVX2FF8MultiplyAdd,
    AVX2FF16Multiply,
    AVX2FF16MultiplyAdd,
    AVX2XorMemory
};

const Ops* InitializeAVX2(const InitializeArgs& args)
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
    return &AVX2Ops;
}

}} // namespace leopard::backend
