/*
    Copyright (c) 2017 Christopher A. Taylor.  All rights reserved.
    See Leopard2Backend.h for the full BSD license notice.
*/

#include "Leopard2Backend.h"

#include <cstring>
#include <memory>
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

static uint64_t ScalarFF8PackedProduct(
    uint64_t value,
    const uint8_t* table)
{
    uint64_t product = 0;
    product ^= static_cast<uint64_t>(table[value & 255U]);
    product ^= static_cast<uint64_t>(table[(value >> 8) & 255U]) << 8;
    product ^= static_cast<uint64_t>(table[(value >> 16) & 255U]) << 16;
    product ^= static_cast<uint64_t>(table[(value >> 24) & 255U]) << 24;
    product ^= static_cast<uint64_t>(table[(value >> 32) & 255U]) << 32;
    product ^= static_cast<uint64_t>(table[(value >> 40) & 255U]) << 40;
    product ^= static_cast<uint64_t>(table[(value >> 48) & 255U]) << 48;
    product ^= static_cast<uint64_t>(table[value >> 56]) << 56;
    return product;
}

template<bool Inverse>
static void ScalarFF8Butterfly2(
    void* x_pointer,
    void* y_pointer,
    uint16_t multiplier_log,
    uint64_t byte_count)
{
    uint8_t* x = static_cast<uint8_t*>(x_pointer);
    uint8_t* y = static_cast<uint8_t*>(y_pointer);
    const uint8_t* table = FF8Table +
        static_cast<size_t>(multiplier_log) * 256U;
    while (byte_count >= 8)
    {
        uint64_t x_value;
        uint64_t y_value;
        std::memcpy(&x_value, x, sizeof(x_value));
        std::memcpy(&y_value, y, sizeof(y_value));
        if (Inverse)
        {
            y_value ^= x_value;
            x_value ^= ScalarFF8PackedProduct(y_value, table);
        }
        else
        {
            x_value ^= ScalarFF8PackedProduct(y_value, table);
            y_value ^= x_value;
        }
        std::memcpy(x, &x_value, sizeof(x_value));
        std::memcpy(y, &y_value, sizeof(y_value));
        x += 8;
        y += 8;
        byte_count -= 8;
    }
    while (byte_count-- != 0)
    {
        if (Inverse)
        {
            *y ^= *x;
            *x ^= table[*y];
        }
        else
        {
            *x ^= table[*y];
            *y ^= *x;
        }
        ++x;
        ++y;
    }
}

static void ScalarFF8IFFTButterfly2(
    void* x, void* y, uint16_t multiplier_log, uint64_t byte_count)
{
    ScalarFF8Butterfly2<true>(x, y, multiplier_log, byte_count);
}

static void ScalarFF8FFTButterfly2(
    void* x, void* y, uint16_t multiplier_log, uint64_t byte_count)
{
    ScalarFF8Butterfly2<false>(x, y, multiplier_log, byte_count);
}

static void ScalarFF8IFFTButterfly2Xor(
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
    const uint8_t* table = FF8Table +
        static_cast<size_t>(multiplier_log) * 256U;
    while (byte_count >= 8)
    {
        uint64_t x_value;
        uint64_t y_value;
        uint64_t x_accumulator;
        uint64_t y_accumulator;
        std::memcpy(&x_value, x_input, sizeof(x_value));
        std::memcpy(&y_value, y_input, sizeof(y_value));
        std::memcpy(&x_accumulator, x_output, sizeof(x_accumulator));
        std::memcpy(&y_accumulator, y_output, sizeof(y_accumulator));
        y_value ^= x_value;
        x_value ^= ScalarFF8PackedProduct(y_value, table);
        x_accumulator ^= x_value;
        y_accumulator ^= y_value;
        std::memcpy(x_output, &x_accumulator, sizeof(x_accumulator));
        std::memcpy(y_output, &y_accumulator, sizeof(y_accumulator));
        x_input += 8;
        y_input += 8;
        x_output += 8;
        y_output += 8;
        byte_count -= 8;
    }
    while (byte_count-- != 0)
    {
        const uint8_t y_value = static_cast<uint8_t>(*y_input ^ *x_input);
        const uint8_t x_value = static_cast<uint8_t>(
            *x_input ^ table[y_value]);
        *x_output ^= x_value;
        *y_output ^= y_value;
        ++x_input;
        ++y_input;
        ++x_output;
        ++y_output;
    }
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

template<bool Inverse>
static void ScalarFF16Butterfly2(
    void* x_pointer,
    void* y_pointer,
    uint16_t multiplier_log,
    uint64_t byte_count)
{
    uint8_t* x = static_cast<uint8_t*>(x_pointer);
    uint8_t* y = static_cast<uint8_t*>(y_pointer);
    uint64_t offset = 0;
    while (byte_count - offset >= 64)
    {
        for (unsigned i = 0; i < 32; ++i)
        {
            uint16_t x_value = static_cast<uint16_t>(x[offset + i] |
                (static_cast<unsigned>(x[offset + 32 + i]) << 8));
            uint16_t y_value = static_cast<uint16_t>(y[offset + i] |
                (static_cast<unsigned>(y[offset + 32 + i]) << 8));
            if (Inverse)
            {
                y_value ^= x_value;
                x_value ^= ScalarFF16Product(multiplier_log, y_value);
            }
            else
            {
                x_value ^= ScalarFF16Product(multiplier_log, y_value);
                y_value ^= x_value;
            }
            x[offset + i] = static_cast<uint8_t>(x_value);
            x[offset + 32 + i] = static_cast<uint8_t>(x_value >> 8);
            y[offset + i] = static_cast<uint8_t>(y_value);
            y[offset + 32 + i] = static_cast<uint8_t>(y_value >> 8);
        }
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
            x_value ^= ScalarFF16Product(multiplier_log, y_value);
        }
        else
        {
            x_value ^= ScalarFF16Product(multiplier_log, y_value);
            y_value ^= x_value;
        }
        x[offset + i] = static_cast<uint8_t>(x_value);
        x[offset + symbols + i] = static_cast<uint8_t>(x_value >> 8);
        y[offset + i] = static_cast<uint8_t>(y_value);
        y[offset + symbols + i] = static_cast<uint8_t>(y_value >> 8);
    }
}

static void ScalarFF16IFFTButterfly2(
    void* x, void* y, uint16_t multiplier_log, uint64_t byte_count)
{
    ScalarFF16Butterfly2<true>(x, y, multiplier_log, byte_count);
}

static void ScalarFF16FFTButterfly2(
    void* x, void* y, uint16_t multiplier_log, uint64_t byte_count)
{
    ScalarFF16Butterfly2<false>(x, y, multiplier_log, byte_count);
}

static void ScalarXorMemory(
    void* destination,
    const void* source,
    uint64_t byte_count)
{
    leopard::xor_mem_baseline(destination, source, byte_count);
}

static void ScalarXorMemory4(
    void* destination0, const void* source0,
    void* destination1, const void* source1,
    void* destination2, const void* source2,
    void* destination3, const void* source3,
    uint64_t byte_count)
{
    ScalarXorMemory(destination0, source0, byte_count);
    ScalarXorMemory(destination1, source1, byte_count);
    ScalarXorMemory(destination2, source2, byte_count);
    ScalarXorMemory(destination3, source3, byte_count);
}

template<bool Inverse>
static void ScalarFF8Butterfly4(
    void* value0, void* value1, void* value2, void* value3,
    uint16_t log01, uint16_t log23, uint16_t log02,
    uint64_t byte_count)
{
    static const uint16_t kZeroSkew = 255;
    if (Inverse)
    {
        if (log01 == kZeroSkew)
            ScalarXorMemory(value1, value0, byte_count);
        else
            ScalarFF8Butterfly2<true>(
                value0, value1, log01, byte_count);
        if (log23 == kZeroSkew)
            ScalarXorMemory(value3, value2, byte_count);
        else
            ScalarFF8Butterfly2<true>(
                value2, value3, log23, byte_count);
        if (log02 == kZeroSkew)
        {
            ScalarXorMemory(value2, value0, byte_count);
            ScalarXorMemory(value3, value1, byte_count);
        }
        else
        {
            ScalarFF8Butterfly2<true>(
                value0, value2, log02, byte_count);
            ScalarFF8Butterfly2<true>(
                value1, value3, log02, byte_count);
        }
    }
    else
    {
        if (log02 == kZeroSkew)
        {
            ScalarXorMemory(value2, value0, byte_count);
            ScalarXorMemory(value3, value1, byte_count);
        }
        else
        {
            ScalarFF8Butterfly2<false>(
                value0, value2, log02, byte_count);
            ScalarFF8Butterfly2<false>(
                value1, value3, log02, byte_count);
        }
        if (log01 == kZeroSkew)
            ScalarXorMemory(value1, value0, byte_count);
        else
            ScalarFF8Butterfly2<false>(
                value0, value1, log01, byte_count);
        if (log23 == kZeroSkew)
            ScalarXorMemory(value3, value2, byte_count);
        else
            ScalarFF8Butterfly2<false>(
                value2, value3, log23, byte_count);
    }
}

static void ScalarFF8IFFTButterfly4(
    void* value0, void* value1, void* value2, void* value3,
    uint16_t log01, uint16_t log23, uint16_t log02,
    uint64_t byte_count)
{
    ScalarFF8Butterfly4<true>(value0, value1, value2, value3,
        log01, log23, log02, byte_count);
}

static void ScalarFF8FFTButterfly4(
    void* value0, void* value1, void* value2, void* value3,
    uint16_t log01, uint16_t log23, uint16_t log02,
    uint64_t byte_count)
{
    ScalarFF8Butterfly4<false>(value0, value1, value2, value3,
        log01, log23, log02, byte_count);
}

template<bool Inverse>
static void ScalarFF16Butterfly4(
    void* value0, void* value1, void* value2, void* value3,
    uint16_t log01, uint16_t log23, uint16_t log02,
    uint64_t byte_count)
{
    static const uint16_t kZeroSkew = 65535;
    if (Inverse)
    {
        if (log01 == kZeroSkew)
            ScalarXorMemory(value1, value0, byte_count);
        else
            ScalarFF16Butterfly2<true>(
                value0, value1, log01, byte_count);
        if (log23 == kZeroSkew)
            ScalarXorMemory(value3, value2, byte_count);
        else
            ScalarFF16Butterfly2<true>(
                value2, value3, log23, byte_count);
        if (log02 == kZeroSkew)
        {
            ScalarXorMemory(value2, value0, byte_count);
            ScalarXorMemory(value3, value1, byte_count);
        }
        else
        {
            ScalarFF16Butterfly2<true>(
                value0, value2, log02, byte_count);
            ScalarFF16Butterfly2<true>(
                value1, value3, log02, byte_count);
        }
    }
    else
    {
        if (log02 == kZeroSkew)
        {
            ScalarXorMemory(value2, value0, byte_count);
            ScalarXorMemory(value3, value1, byte_count);
        }
        else
        {
            ScalarFF16Butterfly2<false>(
                value0, value2, log02, byte_count);
            ScalarFF16Butterfly2<false>(
                value1, value3, log02, byte_count);
        }
        if (log01 == kZeroSkew)
            ScalarXorMemory(value1, value0, byte_count);
        else
            ScalarFF16Butterfly2<false>(
                value0, value1, log01, byte_count);
        if (log23 == kZeroSkew)
            ScalarXorMemory(value3, value2, byte_count);
        else
            ScalarFF16Butterfly2<false>(
                value2, value3, log23, byte_count);
    }
}

static void ScalarFF16IFFTButterfly4(
    void* value0, void* value1, void* value2, void* value3,
    uint16_t log01, uint16_t log23, uint16_t log02,
    uint64_t byte_count)
{
    ScalarFF16Butterfly4<true>(value0, value1, value2, value3,
        log01, log23, log02, byte_count);
}

static void ScalarFF16FFTButterfly4(
    void* value0, void* value1, void* value2, void* value3,
    uint16_t log01, uint16_t log23, uint16_t log02,
    uint64_t byte_count)
{
    ScalarFF16Butterfly4<false>(value0, value1, value2, value3,
        log01, log23, log02, byte_count);
}

static const Ops ScalarOps = {
    LEO2_BACKEND_SCALAR,
    "scalar",
    ScalarFF8Multiply,
    ScalarFF8MultiplyAdd,
    ScalarFF16Multiply,
    ScalarFF16MultiplyAdd,
    ScalarXorMemory,
    ScalarXorMemory4,
    ScalarFF8IFFTButterfly2,
    ScalarFF8FFTButterfly2,
    ScalarFF8IFFTButterfly2Xor,
    ScalarFF8IFFTButterfly4,
    ScalarFF8FFTButterfly4,
    ScalarFF16IFFTButterfly2,
    ScalarFF16FFTButterfly2,
    ScalarFF16IFFTButterfly4,
    ScalarFF16FFTButterfly4
};

const Ops* InitializeScalar(const InitializeArgs& args)
{
    if (!args.ff8_multiply_log || !args.ff16_multiply_log)
        return NULL;

    // Tables are immutable and intentionally retained for process lifetime,
    // but construction must be failure-atomic.  Build the complete pair in
    // local owners and publish neither pointer until both allocations and all
    // arithmetic have completed.
    if (FF8Table || FF16Table)
        return FF8Table && FF16Table ? &ScalarOps : NULL;
#ifdef LEO2_ENABLE_TEST_HOOKS
    if (TestShouldFailAllocation(LEO2_BACKEND_SCALAR, false))
        return NULL;
#endif
    std::unique_ptr<uint8_t[]> ff8(
        new (std::nothrow) uint8_t[256U * 256U]);
    if (!ff8)
        return NULL;
#ifdef LEO2_ENABLE_TEST_HOOKS
    if (TestShouldFailAllocation(LEO2_BACKEND_SCALAR, true))
        return NULL;
#endif
    std::unique_ptr<uint16_t[]> ff16(
        new (std::nothrow) uint16_t[65536U * 64U]);
    if (!ff16)
        return NULL;

    for (unsigned log = 0; log < 256; ++log)
        for (unsigned value = 0; value < 256; ++value)
            ff8[log * 256U + value] = args.ff8_multiply_log(
                static_cast<uint8_t>(value), static_cast<uint8_t>(log));

    for (unsigned log = 0; log < 65536; ++log)
    {
        uint16_t* table = ff16.get() + static_cast<size_t>(log) * 64U;
        for (unsigned nibble = 0; nibble < 4; ++nibble)
            for (unsigned value = 0; value < 16; ++value)
                table[nibble * 16U + value] = args.ff16_multiply_log(
                    static_cast<uint16_t>(value << (nibble * 4)),
                    static_cast<uint16_t>(log));
    }
    FF8Table = ff8.release();
    FF16Table = ff16.release();
    return &ScalarOps;
}

#ifdef LEO2_ENABLE_TEST_HOOKS
void TestGetScalarTableState(TestBackendState* state)
{
    state->ff8_published = FF8Table != NULL;
    state->ff16_published = FF16Table != NULL;
    state->ff8_bytes = 256U * 256U * sizeof(uint8_t);
    state->ff16_bytes = 65536U * 64U * sizeof(uint16_t);
}
#endif

}} // namespace leopard::backend
