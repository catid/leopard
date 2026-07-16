/*
    Copyright (c) 2017 Christopher A. Taylor.  All rights reserved.
    See Leopard2Backend.h for the full BSD license notice.
*/

#include "Leopard2Backend.h"

#include <cstring>
#include <new>

namespace leopard {
void xor_mem_baseline(void* destination, const void* source, uint64_t byte_count);
namespace backend {

static uint8_t* FF8Table = NULL;
static uint16_t* FF16Table = NULL;

template<bool Add>
static void ScalarFF8Operation(
    void* destination,
    const void* source,
    uint16_t multiplier_log,
    uint64_t byte_count)
{
    uint8_t* output = static_cast<uint8_t*>(destination);
    const uint8_t* input = static_cast<const uint8_t*>(source);
    const uint8_t* table = FF8Table +
        static_cast<size_t>(multiplier_log) * 256U;
    while (byte_count >= 8)
    {
        uint64_t result = 0;
        if (Add)
            std::memcpy(&result, output, sizeof(result));
        result ^= static_cast<uint64_t>(table[input[0]]);
        result ^= static_cast<uint64_t>(table[input[1]]) << 8;
        result ^= static_cast<uint64_t>(table[input[2]]) << 16;
        result ^= static_cast<uint64_t>(table[input[3]]) << 24;
        result ^= static_cast<uint64_t>(table[input[4]]) << 32;
        result ^= static_cast<uint64_t>(table[input[5]]) << 40;
        result ^= static_cast<uint64_t>(table[input[6]]) << 48;
        result ^= static_cast<uint64_t>(table[input[7]]) << 56;
        std::memcpy(output, &result, sizeof(result));
        output += 8;
        input += 8;
        byte_count -= 8;
    }
    while (byte_count-- != 0)
    {
        const uint8_t product = table[*input++];
        if (Add)
            *output++ ^= product;
        else
            *output++ = product;
    }
}

static void ScalarFF8Multiply(
    void* destination,
    const void* source,
    uint16_t multiplier_log,
    uint64_t byte_count)
{
    ScalarFF8Operation<false>(
        destination, source, multiplier_log, byte_count);
}

static void ScalarFF8MultiplyAdd(
    void* destination,
    const void* source,
    uint16_t multiplier_log,
    uint64_t byte_count)
{
    ScalarFF8Operation<true>(
        destination, source, multiplier_log, byte_count);
}

static uint16_t ScalarFF16Product(uint16_t multiplier_log, uint16_t value)
{
    const uint16_t* table = FF16Table +
        static_cast<size_t>(multiplier_log) * 64U;
    return static_cast<uint16_t>(
        table[value & 15U] ^
        table[((value >> 4) & 15U) + 16U] ^
        table[((value >> 8) & 15U) + 32U] ^
        table[(value >> 12) + 48U]);
}

template<bool Add>
static void ScalarFF16Operation(
    void* destination,
    const void* source,
    uint16_t multiplier_log,
    uint64_t byte_count)
{
    uint8_t* output = static_cast<uint8_t*>(destination);
    const uint8_t* input = static_cast<const uint8_t*>(source);
    uint64_t offset = 0;
    while (byte_count - offset >= 64)
    {
        for (unsigned i = 0; i < 32; ++i)
        {
            const uint16_t value = static_cast<uint16_t>(input[offset + i] |
                (static_cast<unsigned>(input[offset + 32 + i]) << 8));
            const uint16_t product = ScalarFF16Product(multiplier_log, value);
            if (Add)
            {
                output[offset + i] ^= static_cast<uint8_t>(product);
                output[offset + 32 + i] ^=
                    static_cast<uint8_t>(product >> 8);
            }
            else
            {
                output[offset + i] = static_cast<uint8_t>(product);
                output[offset + 32 + i] = static_cast<uint8_t>(product >> 8);
            }
        }
        offset += 64;
    }

    const uint64_t residual = byte_count - offset;
    const uint64_t symbols = residual / 2;
    for (uint64_t i = 0; i < symbols; ++i)
    {
        const uint16_t value = static_cast<uint16_t>(input[offset + i] |
            (static_cast<unsigned>(input[offset + symbols + i]) << 8));
        const uint16_t product = ScalarFF16Product(multiplier_log, value);
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

static void ScalarFF16Multiply(
    void* destination,
    const void* source,
    uint16_t multiplier_log,
    uint64_t byte_count)
{
    ScalarFF16Operation<false>(
        destination, source, multiplier_log, byte_count);
}

static void ScalarFF16MultiplyAdd(
    void* destination,
    const void* source,
    uint16_t multiplier_log,
    uint64_t byte_count)
{
    ScalarFF16Operation<true>(
        destination, source, multiplier_log, byte_count);
}

static void ScalarXorMemory(
    void* destination,
    const void* source,
    uint64_t byte_count)
{
    leopard::xor_mem_baseline(destination, source, byte_count);
}

static const Ops ScalarOps = {
    LEO2_BACKEND_SCALAR,
    "scalar",
    ScalarFF8Multiply,
    ScalarFF8MultiplyAdd,
    ScalarFF16Multiply,
    ScalarFF16MultiplyAdd,
    ScalarXorMemory
};

const Ops* InitializeScalar(const InitializeArgs& args)
{
    if (!args.ff8_multiply_log || !args.ff16_multiply_log)
        return NULL;
    if (!FF8Table)
        FF8Table = new (std::nothrow) uint8_t[256U * 256U];
    if (!FF16Table)
        FF16Table = new (std::nothrow) uint16_t[65536U * 64U];
    if (!FF8Table || !FF16Table)
        return NULL;

    for (unsigned log = 0; log < 256; ++log)
        for (unsigned value = 0; value < 256; ++value)
            FF8Table[log * 256U + value] = args.ff8_multiply_log(
                static_cast<uint8_t>(value), static_cast<uint8_t>(log));

    for (unsigned log = 0; log < 65536; ++log)
    {
        uint16_t* table = FF16Table + static_cast<size_t>(log) * 64U;
        for (unsigned nibble = 0; nibble < 4; ++nibble)
            for (unsigned value = 0; value < 16; ++value)
                table[nibble * 16U + value] = args.ff16_multiply_log(
                    static_cast<uint16_t>(value << (nibble * 4)),
                    static_cast<uint16_t>(log));
    }
    return &ScalarOps;
}

}} // namespace leopard::backend
