/*
    Copyright (c) 2017 Christopher A. Taylor.  All rights reserved.
    See Leopard2Backend.h for the full BSD license notice.
*/

#include "Leopard2Backend.h"
#include "LeopardCommon.h"

#include <cstring>
#include <memory>
#include <new>

#if defined(__aarch64__)
#include <arm_neon.h>
#endif

namespace leopard {
void xor_mem_baseline(void* destination, const void* source, uint64_t byte_count);
namespace backend {

#ifdef LEO_HAS_FF8
static uint8_t* FF8Table = NULL;
#endif
#ifdef LEO_HAS_FF16
static uint16_t* FF16Table = NULL;
#endif

#ifdef LEO_HAS_FF8
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

static void ScalarFF8FFTButterfly2Out(
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
    const uint8_t* table = multiplier_log == kZeroSkew ? NULL :
        FF8Table + static_cast<size_t>(multiplier_log) * 256U;
    while (byte_count >= 8)
    {
        uint64_t x_value;
        uint64_t y_value;
        std::memcpy(&x_value, x_input, sizeof(x_value));
        std::memcpy(&y_value, y_input, sizeof(y_value));
        if (table)
            x_value ^= ScalarFF8PackedProduct(y_value, table);
        y_value ^= x_value;
        std::memcpy(x_output, &x_value, sizeof(x_value));
        std::memcpy(y_output, &y_value, sizeof(y_value));
        x_input += 8;
        y_input += 8;
        x_output += 8;
        y_output += 8;
        byte_count -= 8;
    }
    while (byte_count-- != 0)
    {
        uint8_t x_value = *x_input++;
        uint8_t y_value = *y_input++;
        if (table)
            x_value ^= table[y_value];
        y_value ^= x_value;
        *x_output++ = x_value;
        *y_output++ = y_value;
    }
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

#endif // LEO_HAS_FF8

#ifdef LEO_HAS_FF16
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

static void ScalarFF16FFTButterfly2Out(
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
    uint64_t offset = 0;
    while (byte_count - offset >= 64)
    {
        for (unsigned i = 0; i < 32; ++i)
        {
            uint16_t x_value = static_cast<uint16_t>(x_input[offset + i] |
                (static_cast<unsigned>(x_input[offset + 32 + i]) << 8));
            uint16_t y_value = static_cast<uint16_t>(y_input[offset + i] |
                (static_cast<unsigned>(y_input[offset + 32 + i]) << 8));
            if (multiplier_log != kZeroSkew)
                x_value ^= ScalarFF16Product(multiplier_log, y_value);
            y_value ^= x_value;
            x_output[offset + i] = static_cast<uint8_t>(x_value);
            x_output[offset + 32 + i] = static_cast<uint8_t>(x_value >> 8);
            y_output[offset + i] = static_cast<uint8_t>(y_value);
            y_output[offset + 32 + i] = static_cast<uint8_t>(y_value >> 8);
        }
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
            x_value ^= ScalarFF16Product(multiplier_log, y_value);
        y_value ^= x_value;
        x_output[offset + i] = static_cast<uint8_t>(x_value);
        x_output[offset + symbols + i] = static_cast<uint8_t>(x_value >> 8);
        y_output[offset + i] = static_cast<uint8_t>(y_value);
        y_output[offset + symbols + i] = static_cast<uint8_t>(y_value >> 8);
    }
}

static void ScalarFF16IFFTButterfly2Xor(
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
    uint64_t offset = 0;
    while (byte_count - offset >= 64)
    {
        for (unsigned i = 0; i < 32; ++i)
        {
            uint16_t x_value = static_cast<uint16_t>(x_input[offset + i] |
                (static_cast<unsigned>(x_input[offset + 32 + i]) << 8));
            uint16_t y_value = static_cast<uint16_t>(y_input[offset + i] |
                (static_cast<unsigned>(y_input[offset + 32 + i]) << 8));
            y_value ^= x_value;
            x_value ^= ScalarFF16Product(multiplier_log, y_value);
            x_output[offset + i] ^= static_cast<uint8_t>(x_value);
            x_output[offset + 32 + i] ^=
                static_cast<uint8_t>(x_value >> 8);
            y_output[offset + i] ^= static_cast<uint8_t>(y_value);
            y_output[offset + 32 + i] ^=
                static_cast<uint8_t>(y_value >> 8);
        }
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
        x_value ^= ScalarFF16Product(multiplier_log, y_value);
        x_output[offset + i] ^= static_cast<uint8_t>(x_value);
        x_output[offset + symbols + i] ^=
            static_cast<uint8_t>(x_value >> 8);
        y_output[offset + i] ^= static_cast<uint8_t>(y_value);
        y_output[offset + symbols + i] ^=
            static_cast<uint8_t>(y_value >> 8);
    }
}

#endif // LEO_HAS_FF16

static void ScalarXorMemory(
    void* destination,
    const void* source,
    uint64_t byte_count)
{
    leopard::xor_mem_baseline(destination, source, byte_count);
}

static void ScalarXorMemory2To1(
    void* destination,
    const void* source0,
    const void* source1,
    uint64_t byte_count)
{
    uint8_t* output = static_cast<uint8_t*>(destination);
    const uint8_t* input0 = static_cast<const uint8_t*>(source0);
    const uint8_t* input1 = static_cast<const uint8_t*>(source1);
#if defined(__aarch64__)
    while (byte_count >= 64)
    {
        for (unsigned offset = 0; offset < 64; offset += 16)
        {
            uint8x16_t result = veorq_u8(
                vld1q_u8(output + offset), vld1q_u8(input0 + offset));
            result = veorq_u8(result, vld1q_u8(input1 + offset));
            vst1q_u8(output + offset, result);
        }
        output += 64;
        input0 += 64;
        input1 += 64;
        byte_count -= 64;
    }
    while (byte_count >= 16)
    {
        uint8x16_t result = veorq_u8(
            vld1q_u8(output), vld1q_u8(input0));
        result = veorq_u8(result, vld1q_u8(input1));
        vst1q_u8(output, result);
        output += 16;
        input0 += 16;
        input1 += 16;
        byte_count -= 16;
    }
#endif
    while (byte_count >= sizeof(uint64_t))
    {
        uint64_t destination_word;
        uint64_t source0_word;
        uint64_t source1_word;
        std::memcpy(&destination_word, output, sizeof(destination_word));
        std::memcpy(&source0_word, input0, sizeof(source0_word));
        std::memcpy(&source1_word, input1, sizeof(source1_word));
        destination_word ^= source0_word ^ source1_word;
        std::memcpy(output, &destination_word, sizeof(destination_word));
        output += sizeof(uint64_t);
        input0 += sizeof(uint64_t);
        input1 += sizeof(uint64_t);
        byte_count -= sizeof(uint64_t);
    }
    while (byte_count-- != 0)
        *output++ ^= *input0++ ^ *input1++;
}

static void ScalarXorMemorySources(
    void* destination,
    const void* initial_source,
    const void* const* sources,
    uint32_t source_count,
    uint64_t byte_count)
{
    std::memcpy(destination, initial_source, static_cast<size_t>(byte_count));
    const void* waiting[4];
    unsigned waiting_count = 0;
    for (uint32_t i = 0; i < source_count; ++i)
    {
        if (!sources[i])
            continue;
        waiting[waiting_count++] = sources[i];
        if (waiting_count == 4)
        {
            uint8_t* output = static_cast<uint8_t*>(destination);
            const uint8_t* const input0 =
                static_cast<const uint8_t*>(waiting[0]);
            const uint8_t* const input1 =
                static_cast<const uint8_t*>(waiting[1]);
            const uint8_t* const input2 =
                static_cast<const uint8_t*>(waiting[2]);
            const uint8_t* const input3 =
                static_cast<const uint8_t*>(waiting[3]);
            uint64_t offset = 0;
            while (byte_count - offset >= sizeof(uint64_t))
            {
                uint64_t result;
                uint64_t value0;
                uint64_t value1;
                uint64_t value2;
                uint64_t value3;
                std::memcpy(&result, output + offset, sizeof(result));
                std::memcpy(&value0, input0 + offset, sizeof(value0));
                std::memcpy(&value1, input1 + offset, sizeof(value1));
                std::memcpy(&value2, input2 + offset, sizeof(value2));
                std::memcpy(&value3, input3 + offset, sizeof(value3));
                result ^= value0 ^ value1 ^ value2 ^ value3;
                std::memcpy(output + offset, &result, sizeof(result));
                offset += sizeof(uint64_t);
            }
            while (offset < byte_count)
            {
                output[offset] ^= input0[offset] ^ input1[offset] ^
                    input2[offset] ^ input3[offset];
                ++offset;
            }
            waiting_count = 0;
        }
    }
    if (waiting_count >= 2)
    {
        ScalarXorMemory2To1(
            destination, waiting[0], waiting[1], byte_count);
        waiting_count -= 2;
        if (waiting_count != 0)
            waiting[0] = waiting[2];
    }
    if (waiting_count != 0)
        ScalarXorMemory(destination, waiting[0], byte_count);
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

#ifdef LEO_HAS_FF8
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

static void ScalarFF8IFFTButterfly4Out(
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
    const uint8_t* table01 = log01 == kZeroSkew ? NULL :
        FF8Table + static_cast<size_t>(log01) * 256U;
    const uint8_t* table23 = log23 == kZeroSkew ? NULL :
        FF8Table + static_cast<size_t>(log23) * 256U;
    const uint8_t* table02 = log02 == kZeroSkew ? NULL :
        FF8Table + static_cast<size_t>(log02) * 256U;
    while (byte_count >= 8)
    {
        uint64_t x0, x1, x2, x3;
        std::memcpy(&x0, input0, sizeof(x0));
        std::memcpy(&x1, input1, sizeof(x1));
        std::memcpy(&x2, input2, sizeof(x2));
        std::memcpy(&x3, input3, sizeof(x3));
        x1 ^= x0;
        if (table01)
            x0 ^= ScalarFF8PackedProduct(x1, table01);
        x3 ^= x2;
        if (table23)
            x2 ^= ScalarFF8PackedProduct(x3, table23);
        x2 ^= x0;
        x3 ^= x1;
        if (table02)
        {
            x0 ^= ScalarFF8PackedProduct(x2, table02);
            x1 ^= ScalarFF8PackedProduct(x3, table02);
        }
        std::memcpy(output0, &x0, sizeof(x0));
        std::memcpy(output1, &x1, sizeof(x1));
        std::memcpy(output2, &x2, sizeof(x2));
        std::memcpy(output3, &x3, sizeof(x3));
        input0 += 8;
        input1 += 8;
        input2 += 8;
        input3 += 8;
        output0 += 8;
        output1 += 8;
        output2 += 8;
        output3 += 8;
        byte_count -= 8;
    }
    while (byte_count-- != 0)
    {
        uint8_t x0 = *input0++;
        uint8_t x1 = *input1++;
        uint8_t x2 = *input2++;
        uint8_t x3 = *input3++;
        x1 ^= x0;
        if (table01)
            x0 ^= table01[x1];
        x3 ^= x2;
        if (table23)
            x2 ^= table23[x3];
        x2 ^= x0;
        x3 ^= x1;
        if (table02)
        {
            x0 ^= table02[x2];
            x1 ^= table02[x3];
        }
        *output0++ = x0;
        *output1++ = x1;
        *output2++ = x2;
        *output3++ = x3;
    }
}

static void ScalarFF8FFTButterfly4Out(
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
    const uint8_t* table01 = log01 == kZeroSkew ? NULL :
        FF8Table + static_cast<size_t>(log01) * 256U;
    const uint8_t* table23 = log23 == kZeroSkew ? NULL :
        FF8Table + static_cast<size_t>(log23) * 256U;
    const uint8_t* table02 = log02 == kZeroSkew ? NULL :
        FF8Table + static_cast<size_t>(log02) * 256U;
    while (byte_count >= 8)
    {
        uint64_t x0, x1, x2, x3;
        std::memcpy(&x0, input0, sizeof(x0));
        std::memcpy(&x1, input1, sizeof(x1));
        std::memcpy(&x2, input2, sizeof(x2));
        std::memcpy(&x3, input3, sizeof(x3));
        if (table02)
        {
            x0 ^= ScalarFF8PackedProduct(x2, table02);
            x1 ^= ScalarFF8PackedProduct(x3, table02);
        }
        x2 ^= x0;
        x3 ^= x1;
        if (table01)
            x0 ^= ScalarFF8PackedProduct(x1, table01);
        x1 ^= x0;
        if (table23)
            x2 ^= ScalarFF8PackedProduct(x3, table23);
        x3 ^= x2;
        std::memcpy(output0, &x0, sizeof(x0));
        std::memcpy(output1, &x1, sizeof(x1));
        std::memcpy(output2, &x2, sizeof(x2));
        std::memcpy(output3, &x3, sizeof(x3));
        input0 += 8;
        input1 += 8;
        input2 += 8;
        input3 += 8;
        output0 += 8;
        output1 += 8;
        output2 += 8;
        output3 += 8;
        byte_count -= 8;
    }
    while (byte_count-- != 0)
    {
        uint8_t x0 = *input0++;
        uint8_t x1 = *input1++;
        uint8_t x2 = *input2++;
        uint8_t x3 = *input3++;
        if (table02)
        {
            x0 ^= table02[x2];
            x1 ^= table02[x3];
        }
        x2 ^= x0;
        x3 ^= x1;
        if (table01)
            x0 ^= table01[x1];
        x1 ^= x0;
        if (table23)
            x2 ^= table23[x3];
        x3 ^= x2;
        *output0++ = x0;
        *output1++ = x1;
        *output2++ = x2;
        *output3++ = x3;
    }
}

template<bool Inverse>
static void ScalarFF8Butterfly4Range(
    void* const* work, unsigned distance,
    uint16_t log01, uint16_t log23, uint16_t log02,
    uint64_t byte_count, bool prefer_fused)
{
    (void)prefer_fused;
    for (unsigned i = 0; i < distance; ++i)
    {
        ScalarFF8Butterfly4<Inverse>(
            work[i], work[i + distance],
            work[i + distance * 2U], work[i + distance * 3U],
            log01, log23, log02, byte_count);
    }
}

static void ScalarFF8IFFTButterfly4Range(
    void* const* work, unsigned distance,
    uint16_t log01, uint16_t log23, uint16_t log02,
    uint64_t byte_count, bool prefer_fused)
{
    ScalarFF8Butterfly4Range<true>(work, distance,
        log01, log23, log02, byte_count, prefer_fused);
}

static void ScalarFF8FFTButterfly4Range(
    void* const* work, unsigned distance,
    uint16_t log01, uint16_t log23, uint16_t log02,
    uint64_t byte_count, bool prefer_fused)
{
    ScalarFF8Butterfly4Range<false>(work, distance,
        log01, log23, log02, byte_count, prefer_fused);
}

static void ScalarFF8IFFTButterfly4Xor(
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
    const uint8_t* table01 = log01 == kZeroSkew ? NULL :
        FF8Table + static_cast<size_t>(log01) * 256U;
    const uint8_t* table23 = log23 == kZeroSkew ? NULL :
        FF8Table + static_cast<size_t>(log23) * 256U;
    const uint8_t* table02 = log02 == kZeroSkew ? NULL :
        FF8Table + static_cast<size_t>(log02) * 256U;
    while (byte_count >= 8)
    {
        uint64_t x0, x1, x2, x3;
        uint64_t accumulator0, accumulator1, accumulator2, accumulator3;
        std::memcpy(&x0, input0, sizeof(x0));
        std::memcpy(&x1, input1, sizeof(x1));
        std::memcpy(&x2, input2, sizeof(x2));
        std::memcpy(&x3, input3, sizeof(x3));
        std::memcpy(&accumulator0, output0, sizeof(accumulator0));
        std::memcpy(&accumulator1, output1, sizeof(accumulator1));
        std::memcpy(&accumulator2, output2, sizeof(accumulator2));
        std::memcpy(&accumulator3, output3, sizeof(accumulator3));
        x1 ^= x0;
        if (table01)
            x0 ^= ScalarFF8PackedProduct(x1, table01);
        x3 ^= x2;
        if (table23)
            x2 ^= ScalarFF8PackedProduct(x3, table23);
        x2 ^= x0;
        x3 ^= x1;
        if (table02)
        {
            x0 ^= ScalarFF8PackedProduct(x2, table02);
            x1 ^= ScalarFF8PackedProduct(x3, table02);
        }
        accumulator0 ^= x0;
        accumulator1 ^= x1;
        accumulator2 ^= x2;
        accumulator3 ^= x3;
        std::memcpy(output0, &accumulator0, sizeof(accumulator0));
        std::memcpy(output1, &accumulator1, sizeof(accumulator1));
        std::memcpy(output2, &accumulator2, sizeof(accumulator2));
        std::memcpy(output3, &accumulator3, sizeof(accumulator3));
        input0 += 8;
        input1 += 8;
        input2 += 8;
        input3 += 8;
        output0 += 8;
        output1 += 8;
        output2 += 8;
        output3 += 8;
        byte_count -= 8;
    }
    while (byte_count-- != 0)
    {
        uint8_t x0 = *input0++;
        uint8_t x1 = *input1++;
        uint8_t x2 = *input2++;
        uint8_t x3 = *input3++;
        x1 ^= x0;
        if (table01)
            x0 ^= table01[x1];
        x3 ^= x2;
        if (table23)
            x2 ^= table23[x3];
        x2 ^= x0;
        x3 ^= x1;
        if (table02)
        {
            x0 ^= table02[x2];
            x1 ^= table02[x3];
        }
        *output0++ ^= x0;
        *output1++ ^= x1;
        *output2++ ^= x2;
        *output3++ ^= x3;
    }
}

static void ScalarFF8IFFTButterfly4XorRange(
    void* const* work, void* const* xor_output, unsigned distance,
    uint16_t log01, uint16_t log23, uint16_t log02,
    uint64_t byte_count)
{
    for (unsigned i = 0; i < distance; ++i)
        ScalarFF8IFFTButterfly4Xor(
            work[i], work[i + distance],
            work[i + distance * 2U], work[i + distance * 3U],
            xor_output[i], xor_output[i + distance],
            xor_output[i + distance * 2U],
            xor_output[i + distance * 3U],
            log01, log23, log02, byte_count);
}

#endif // LEO_HAS_FF8

#ifdef LEO_HAS_FF16
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

template<bool Inverse>
static void ScalarFF16Butterfly4Range(
    void* const* work, unsigned distance,
    uint16_t log01, uint16_t log23, uint16_t log02,
    uint64_t byte_count, bool prefer_fused)
{
    (void)prefer_fused;
    for (unsigned i = 0; i < distance; ++i)
    {
        ScalarFF16Butterfly4<Inverse>(
            work[i], work[i + distance],
            work[i + distance * 2U], work[i + distance * 3U],
            log01, log23, log02, byte_count);
    }
}

static void ScalarFF16IFFTButterfly4Range(
    void* const* work, unsigned distance,
    uint16_t log01, uint16_t log23, uint16_t log02,
    uint64_t byte_count, bool prefer_fused)
{
    ScalarFF16Butterfly4Range<true>(work, distance,
        log01, log23, log02, byte_count, prefer_fused);
}

static void ScalarFF16FFTButterfly4Range(
    void* const* work, unsigned distance,
    uint16_t log01, uint16_t log23, uint16_t log02,
    uint64_t byte_count, bool prefer_fused)
{
    ScalarFF16Butterfly4Range<false>(work, distance,
        log01, log23, log02, byte_count, prefer_fused);
}

static inline void ScalarFF16FFTButterfly4Values(
    uint16_t& x0, uint16_t& x1, uint16_t& x2, uint16_t& x3,
    uint16_t log01, uint16_t log23, uint16_t log02)
{
    static const uint16_t kZeroSkew = 65535;
    if (log02 != kZeroSkew)
    {
        x0 ^= ScalarFF16Product(log02, x2);
        x1 ^= ScalarFF16Product(log02, x3);
    }
    x2 ^= x0;
    x3 ^= x1;
    if (log01 != kZeroSkew)
        x0 ^= ScalarFF16Product(log01, x1);
    x1 ^= x0;
    if (log23 != kZeroSkew)
        x2 ^= ScalarFF16Product(log23, x3);
    x3 ^= x2;
}

static inline void ScalarFF16IFFTButterfly4Values(
    uint16_t& x0, uint16_t& x1, uint16_t& x2, uint16_t& x3,
    uint16_t log01, uint16_t log23, uint16_t log02)
{
    static const uint16_t kZeroSkew = 65535;
    x1 ^= x0;
    if (log01 != kZeroSkew)
        x0 ^= ScalarFF16Product(log01, x1);
    x3 ^= x2;
    if (log23 != kZeroSkew)
        x2 ^= ScalarFF16Product(log23, x3);
    x2 ^= x0;
    x3 ^= x1;
    if (log02 != kZeroSkew)
    {
        x0 ^= ScalarFF16Product(log02, x2);
        x1 ^= ScalarFF16Product(log02, x3);
    }
}

static void ScalarFF16IFFTButterfly4Out(
    const void* input0_pointer, const void* input1_pointer,
    const void* input2_pointer, const void* input3_pointer,
    void* output0_pointer, void* output1_pointer,
    void* output2_pointer, void* output3_pointer,
    uint16_t log01, uint16_t log23, uint16_t log02,
    uint64_t byte_count)
{
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
        for (unsigned i = 0; i < 32; ++i)
        {
            uint16_t values[4];
            for (unsigned lane = 0; lane < 4; ++lane)
                values[lane] = static_cast<uint16_t>(inputs[lane][offset + i] |
                    (static_cast<unsigned>(inputs[lane][offset + 32 + i]) << 8));
            ScalarFF16IFFTButterfly4Values(
                values[0], values[1], values[2], values[3],
                log01, log23, log02);
            for (unsigned lane = 0; lane < 4; ++lane)
            {
                outputs[lane][offset + i] = static_cast<uint8_t>(values[lane]);
                outputs[lane][offset + 32 + i] =
                    static_cast<uint8_t>(values[lane] >> 8);
            }
        }
        offset += 64;
    }
    const uint64_t symbols = (byte_count - offset) / 2;
    for (uint64_t i = 0; i < symbols; ++i)
    {
        uint16_t values[4];
        for (unsigned lane = 0; lane < 4; ++lane)
            values[lane] = static_cast<uint16_t>(inputs[lane][offset + i] |
                (static_cast<unsigned>(inputs[lane][offset + symbols + i]) << 8));
        ScalarFF16IFFTButterfly4Values(
            values[0], values[1], values[2], values[3],
            log01, log23, log02);
        for (unsigned lane = 0; lane < 4; ++lane)
        {
            outputs[lane][offset + i] = static_cast<uint8_t>(values[lane]);
            outputs[lane][offset + symbols + i] =
                static_cast<uint8_t>(values[lane] >> 8);
        }
    }
}

static void ScalarFF16FFTButterfly4Out(
    const void* input0_pointer, const void* input1_pointer,
    const void* input2_pointer, const void* input3_pointer,
    void* output0_pointer, void* output1_pointer,
    void* output2_pointer, void* output3_pointer,
    uint16_t log01, uint16_t log23, uint16_t log02,
    uint64_t byte_count)
{
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
        for (unsigned i = 0; i < 32; ++i)
        {
            uint16_t values[4];
            for (unsigned lane = 0; lane < 4; ++lane)
                values[lane] = static_cast<uint16_t>(inputs[lane][offset + i] |
                    (static_cast<unsigned>(inputs[lane][offset + 32 + i]) << 8));
            ScalarFF16FFTButterfly4Values(
                values[0], values[1], values[2], values[3],
                log01, log23, log02);
            for (unsigned lane = 0; lane < 4; ++lane)
            {
                outputs[lane][offset + i] = static_cast<uint8_t>(values[lane]);
                outputs[lane][offset + 32 + i] =
                    static_cast<uint8_t>(values[lane] >> 8);
            }
        }
        offset += 64;
    }
    const uint64_t symbols = (byte_count - offset) / 2;
    for (uint64_t i = 0; i < symbols; ++i)
    {
        uint16_t values[4];
        for (unsigned lane = 0; lane < 4; ++lane)
            values[lane] = static_cast<uint16_t>(inputs[lane][offset + i] |
                (static_cast<unsigned>(inputs[lane][offset + symbols + i]) << 8));
        ScalarFF16FFTButterfly4Values(
            values[0], values[1], values[2], values[3],
            log01, log23, log02);
        for (unsigned lane = 0; lane < 4; ++lane)
        {
            outputs[lane][offset + i] = static_cast<uint8_t>(values[lane]);
            outputs[lane][offset + symbols + i] =
                static_cast<uint8_t>(values[lane] >> 8);
        }
    }
}
#endif // LEO_HAS_FF16

static const Ops ScalarOps = {
    LEO2_BACKEND_SCALAR,
    "scalar",
#ifdef LEO_HAS_FF8
    ScalarFF8Multiply,
    ScalarFF8MultiplyAdd,
#else
    NULL,
    NULL,
#endif
#ifdef LEO_HAS_FF16
    ScalarFF16Multiply,
    ScalarFF16MultiplyAdd,
#else
    NULL,
    NULL,
#endif
    ScalarXorMemory,
    ScalarXorMemory2To1,
    ScalarXorMemorySources,
    ScalarXorMemory4,
#ifdef LEO_HAS_FF8
    ScalarFF8IFFTButterfly2,
    ScalarFF8FFTButterfly2,
    ScalarFF8FFTButterfly2Out,
    ScalarFF8IFFTButterfly2Xor,
    ScalarFF8IFFTButterfly4,
    ScalarFF8FFTButterfly4,
    ScalarFF8IFFTButterfly4Out,
    NULL,
    ScalarFF8FFTButterfly4Out,
    ScalarFF8IFFTButterfly4Range,
    ScalarFF8FFTButterfly4Range,
    ScalarFF8IFFTButterfly4XorRange,
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
    ScalarFF16IFFTButterfly2,
    ScalarFF16FFTButterfly2,
    ScalarFF16FFTButterfly2Out,
    ScalarFF16IFFTButterfly2Xor,
    ScalarFF16IFFTButterfly4,
    ScalarFF16FFTButterfly4,
    ScalarFF16IFFTButterfly4Out,
    ScalarFF16FFTButterfly4Out,
    ScalarFF16IFFTButterfly4Range,
    ScalarFF16FFTButterfly4Range
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
    , NULL
    , 0
    , NULL
    , NULL
    , NULL
    , NULL
    , NULL
    , NULL
    // ff8_linear_combination2
    , NULL
    , NULL
    // ff8_ifft_butterfly8_out / ff8_fft_butterfly8_out: this backend keeps the
    // radix-four staging schedule and publishes no radix-eight round.
    , NULL
    , NULL
    , NULL
    , NULL
    // copy_memory
    , NULL
    // xor_memory_sources_group4
    , NULL
};

const Ops* InitializeScalar(const InitializeArgs& args)
{
#ifdef LEO_HAS_FF8
    if (!args.ff8_multiply_log)
        return NULL;
#endif
#ifdef LEO_HAS_FF16
    if (!args.ff16_multiply_log)
        return NULL;
#endif

    // Tables are immutable and intentionally retained for process lifetime,
    // but construction must be failure-atomic.  Build every enabled table in
    // local owners and publish none until all allocations and arithmetic have
    // completed.
#if defined(LEO_HAS_FF8) && defined(LEO_HAS_FF16)
    if (FF8Table || FF16Table)
        return FF8Table && FF16Table ? &ScalarOps : NULL;
#elif defined(LEO_HAS_FF8)
    if (FF8Table)
        return &ScalarOps;
#else
    if (FF16Table)
        return &ScalarOps;
#endif
#ifdef LEO_HAS_FF8
#ifdef LEO2_ENABLE_TEST_HOOKS
    if (TestShouldFailAllocation(LEO2_BACKEND_SCALAR, false))
        return NULL;
#endif
    std::unique_ptr<uint8_t[]> ff8(
        new (std::nothrow) uint8_t[256U * 256U]);
    if (!ff8)
        return NULL;
#endif
#ifdef LEO_HAS_FF16
#ifdef LEO2_ENABLE_TEST_HOOKS
    if (TestShouldFailAllocation(LEO2_BACKEND_SCALAR, true))
        return NULL;
#endif
    std::unique_ptr<uint16_t[]> ff16(
        new (std::nothrow) uint16_t[65536U * 64U]);
    if (!ff16)
        return NULL;
#endif

#ifdef LEO_HAS_FF8
    for (unsigned log = 0; log < 256; ++log)
        for (unsigned value = 0; value < 256; ++value)
            ff8[log * 256U + value] = args.ff8_multiply_log(
                static_cast<uint8_t>(value), static_cast<uint8_t>(log));
#endif

#ifdef LEO_HAS_FF16
    for (unsigned log = 0; log < 65536; ++log)
    {
        uint16_t* table = ff16.get() + static_cast<size_t>(log) * 64U;
        for (unsigned nibble = 0; nibble < 4; ++nibble)
            for (unsigned value = 0; value < 16; ++value)
                table[nibble * 16U + value] = args.ff16_multiply_log(
                    static_cast<uint16_t>(value << (nibble * 4)),
                    static_cast<uint16_t>(log));
    }
#endif
#ifdef LEO_HAS_FF8
    FF8Table = ff8.release();
#endif
#ifdef LEO_HAS_FF16
    FF16Table = ff16.release();
#endif
    return &ScalarOps;
}

#ifdef LEO2_ENABLE_TEST_HOOKS
void TestGetScalarTableState(TestBackendState* state)
{
#ifdef LEO_HAS_FF8
    state->ff8_published = FF8Table != NULL;
    state->ff8_bytes = 256U * 256U * sizeof(uint8_t);
#else
    state->ff8_published = false;
    state->ff8_bytes = 0;
#endif
#ifdef LEO_HAS_FF16
    state->ff16_published = FF16Table != NULL;
    state->ff16_bytes = 65536U * 64U * sizeof(uint16_t);
#else
    state->ff16_published = false;
    state->ff16_bytes = 0;
#endif
}
#endif

}} // namespace leopard::backend
