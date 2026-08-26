/*
    Correctness gate for the isolated AVX-512/GFNI K=R=T=16, B=64 leaf.
    The test translation unit remains baseline-safe and calls the raised-ISA
    object only after the ordinary CPUID/XCR0 classifier succeeds.
*/

#include "Leopard2Backend.h"
#include "LeopardFF8.h"
#include "leopard.h"

#include <array>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <stdexcept>

namespace {

static const unsigned kSide = 16;
static const unsigned kShardBytes = 64;
static const size_t kPayloadBytes = kSide * kShardBytes;
static const uint8_t kZeroSkewLog = 255;

void require(bool condition, const char* message)
{
    if (!condition)
        throw std::runtime_error(message);
}

void scalar_transform(
    const uint8_t input[kPayloadBytes],
    uint8_t output[kPayloadBytes])
{
    uint8_t state[kSide][kShardBytes];
    std::memcpy(state, input, sizeof(state));
    const uint8_t* const skew = leopard::ff8::SkewLogTable();
    require(skew != NULL, "GF8 skew table is unavailable");

    for (unsigned level = 0; level < 4; ++level)
    {
        const unsigned distance = 1U << level;
        const unsigned width = distance * 2U;
        for (unsigned row = 0; row < kSide; row += width)
        {
            const uint8_t log = skew[kSide + row + distance];
            for (unsigned lane = 0; lane < distance; ++lane)
            {
                uint8_t* const x = state[row + lane];
                uint8_t* const y = state[row + lane + distance];
                for (unsigned column = 0; column < kShardBytes; ++column)
                {
                    y[column] ^= x[column];
                    if (log != kZeroSkewLog)
                    {
                        x[column] ^= leopard::ff8::MultiplyLogElement(
                            y[column], log);
                    }
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
            const uint8_t log = skew[row + distance];
            for (unsigned lane = 0; lane < distance; ++lane)
            {
                uint8_t* const x = state[row + lane];
                uint8_t* const y = state[row + lane + distance];
                for (unsigned column = 0; column < kShardBytes; ++column)
                {
                    if (log != kZeroSkewLog)
                    {
                        x[column] ^= leopard::ff8::MultiplyLogElement(
                            y[column], log);
                    }
                    y[column] ^= x[column];
                }
            }
        }
    }

    std::memcpy(output, state, sizeof(state));
}

void check_case(
    leopard::backend::FF8K16R16B64PackedKernel kernel,
    const uint8_t input[kPayloadBytes])
{
    static const size_t kGuardBytes = 79;
    std::array<uint8_t, kPayloadBytes> expected;
    std::array<uint8_t, kPayloadBytes> immutable_input;
    std::array<uint8_t, kPayloadBytes + 2 * kGuardBytes + 8> input_storage;
    std::array<uint8_t, kPayloadBytes + 2 * kGuardBytes + 8> output_storage;
    uint8_t* const actual_input = input_storage.data() + kGuardBytes + 1;
    uint8_t* const actual_output = output_storage.data() + kGuardBytes + 3;

    scalar_transform(input, expected.data());
    std::memcpy(immutable_input.data(), input, kPayloadBytes);
    input_storage.fill(0xa7);
    output_storage.fill(0x5c);
    std::memcpy(actual_input, input, kPayloadBytes);
    kernel(actual_input, actual_output);

    require(std::memcmp(actual_output, expected.data(), kPayloadBytes) == 0,
        "T16 candidate disagrees with scalar transform");
    require(std::memcmp(actual_input, immutable_input.data(), kPayloadBytes) == 0,
        "T16 candidate modified its input");
    for (size_t i = 0; i < kGuardBytes + 3; ++i)
        require(output_storage[i] == 0x5c, "T16 candidate underwrote output");
    for (size_t i = kGuardBytes + 3 + kPayloadBytes;
        i < output_storage.size(); ++i)
        require(output_storage[i] == 0x5c, "T16 candidate overwrote output");
    for (size_t i = 0; i < kGuardBytes + 1; ++i)
        require(input_storage[i] == 0xa7, "T16 candidate under-read guard changed");
    for (size_t i = kGuardBytes + 1 + kPayloadBytes;
        i < input_storage.size(); ++i)
        require(input_storage[i] == 0xa7, "T16 candidate input guard changed");
}

void exhaustive_coordinates(
    leopard::backend::FF8K16R16B64PackedKernel kernel)
{
    std::array<uint8_t, kPayloadBytes> input;
    for (unsigned row = 0; row < kSide; ++row)
    {
        for (unsigned value = 0; value < 256; ++value)
        {
            input.fill(0);
            for (unsigned column = 0; column < kShardBytes; ++column)
            {
                input[row * kShardBytes + column] = static_cast<uint8_t>(
                    value + column * 29U);
            }
            check_case(kernel, input.data());
        }
    }
}

void dense_cases(leopard::backend::FF8K16R16B64PackedKernel kernel)
{
    std::array<uint8_t, kPayloadBytes> input;
    uint64_t random = UINT64_C(0x7431362d64656e73);
    for (unsigned trial = 0; trial < 128; ++trial)
    {
        for (size_t i = 0; i < input.size(); ++i)
        {
            random += UINT64_C(0x9e3779b97f4a7c15);
            uint64_t value = random;
            value = (value ^ (value >> 30)) *
                UINT64_C(0xbf58476d1ce4e5b9);
            value = (value ^ (value >> 27)) *
                UINT64_C(0x94d049bb133111eb);
            value ^= value >> 31;
            input[i] = static_cast<uint8_t>(value + trial + i * 17U);
        }
        check_case(kernel, input.data());
    }
}

} // namespace

int main()
{
    try
    {
        require(leo_init() == Leopard_Success, "Leopard initialization failed");
        const leopard::backend::X86Features features =
            leopard::backend::DetectX86Features();
        if (!features.avx512 || !features.gfni)
        {
            std::puts("AVX-512/GFNI T16 prototype skipped: ISA unavailable");
            return 0;
        }

        std::array<uint8_t, 256> corrupted_skew;
        std::memcpy(corrupted_skew.data(), leopard::ff8::SkewLogTable(),
            corrupted_skew.size());
        corrupted_skew[kSide + 1U] ^= 1U;
        require(leopard::backend::InitializeAVX512GFNIT16(
                    leopard::ff8::MultiplyLogElement,
                    corrupted_skew.data()) == NULL,
            "AVX-512/GFNI T16 KAT accepted a corrupted skew table");

        const leopard::backend::FF8K16R16B64PackedKernel kernel =
            leopard::backend::InitializeAVX512GFNIT16(
                leopard::ff8::MultiplyLogElement,
                leopard::ff8::SkewLogTable());
        require(kernel != NULL, "AVX-512/GFNI T16 private KAT failed");
        exhaustive_coordinates(kernel);
        dense_cases(kernel);
        std::puts("AVX-512/GFNI T16 prototype correctness passed");
        return 0;
    }
    catch (const std::exception& error)
    {
        std::fprintf(stderr, "%s\n", error.what());
        return 1;
    }
}
