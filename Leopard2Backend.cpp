/*
    Copyright (c) 2017 Christopher A. Taylor.  All rights reserved.
    See Leopard2Backend.h for the full BSD license notice.
*/

#include "Leopard2Backend.h"
#include "LeopardCommon.h"

#include <atomic>
#include <cstring>
#include <mutex>

namespace leopard { namespace backend {

static const Ops* SelectedOps = NULL;
static const Ops* QualifiedOps[LEO2_BACKEND_NEON + 1] = {};
enum QualificationState
{
    QualificationUnattempted,
    QualificationPassed,
    QualificationFailed
};
static QualificationState QualificationStates[LEO2_BACKEND_NEON + 1] = {};
static QualificationStatus QualificationFailures[LEO2_BACKEND_NEON + 1] = {};
static InitializeArgs SavedInitializeArgs = { NULL, NULL };
static uint32_t QualifiableBackendMask = 0;
static bool SelfTestPassed = false;
static QualificationStatus StartupFailure = QualificationAvailable;

#ifdef LEO2_ENABLE_TEST_HOOKS
static std::atomic<unsigned> TestFault(TestSetupFaultNone);
static std::atomic<unsigned> TestFaultConsumptions(0);

static TestSetupFault AllocationFaultFor(
    leo2_backend backend,
    bool ff16)
{
    switch (backend)
    {
    case LEO2_BACKEND_SCALAR:
        return ff16 ? TestSetupFaultScalarFF16Allocation
                    : TestSetupFaultScalarFF8Allocation;
    case LEO2_BACKEND_SSSE3:
        return ff16 ? TestSetupFaultSSSE3FF16Allocation
                    : TestSetupFaultSSSE3FF8Allocation;
    case LEO2_BACKEND_AVX2:
        return ff16 ? TestSetupFaultAVX2FF16Allocation
                    : TestSetupFaultAVX2FF8Allocation;
    default:
        return TestSetupFaultNone;
    }
}

static TestSetupFault KATFaultFor(leo2_backend backend)
{
    switch (backend)
    {
    case LEO2_BACKEND_SCALAR: return TestSetupFaultScalarKAT;
    case LEO2_BACKEND_SSSE3: return TestSetupFaultSSSE3KAT;
    case LEO2_BACKEND_AVX2: return TestSetupFaultAVX2KAT;
    default: return TestSetupFaultNone;
    }
}

static bool ConsumeTestFault(TestSetupFault expected)
{
    if (expected == TestSetupFaultNone)
        return false;
    unsigned value = static_cast<unsigned>(expected);
    if (!TestFault.compare_exchange_strong(value,
            static_cast<unsigned>(TestSetupFaultNone),
            std::memory_order_acq_rel, std::memory_order_acquire))
        return false;
    TestFaultConsumptions.fetch_add(1, std::memory_order_relaxed);
    return true;
}
#endif

#ifdef LEO_HAS_FF8
static bool TestFF8(const Ops& ops, FF8MultiplyLog reference)
{
    static const size_t kBytes = 521;
    uint8_t source[kBytes + 2];
    uint8_t output[kBytes + 2];
    uint8_t expected[kBytes + 2];
    for (size_t i = 0; i < kBytes + 2; ++i)
        source[i] = static_cast<uint8_t>((i * 73U + 19U) & 255U);

    for (unsigned log = 0; log < 256; ++log)
    {
        std::memset(output, 0xa5, sizeof(output));
        std::memset(expected, 0xa5, sizeof(expected));
        for (size_t i = 0; i < kBytes; ++i)
            expected[i + 1] = reference(
                source[i + 1], static_cast<uint8_t>(log));
        ops.ff8_multiply(
            output + 1, source + 1, static_cast<uint16_t>(log), kBytes);
        if (std::memcmp(output, expected, sizeof(output)) != 0)
            return false;

        for (size_t i = 0; i < sizeof(output); ++i)
            output[i] = expected[i] =
                static_cast<uint8_t>((i * 29U + log) & 255U);
        for (size_t i = 0; i < kBytes; ++i)
            expected[i + 1] ^= reference(
                source[i + 1], static_cast<uint8_t>(log));
        ops.ff8_multiply_add(
            output + 1, source + 1, static_cast<uint16_t>(log), kBytes);
        if (std::memcmp(output, expected, sizeof(output)) != 0)
            return false;
    }
    return true;
}

template<bool Inverse>
static void ReferenceFF8Butterfly2(
    uint8_t* x,
    uint8_t* y,
    uint8_t log,
    uint64_t byte_count,
    FF8MultiplyLog reference)
{
    for (uint64_t i = 0; i < byte_count; ++i)
    {
        if (Inverse)
        {
            y[i] ^= x[i];
            x[i] ^= reference(y[i], log);
        }
        else
        {
            x[i] ^= reference(y[i], log);
            y[i] ^= x[i];
        }
    }
}

static bool TestFF8Butterflies(const Ops& ops, FF8MultiplyLog reference)
{
    static const uint64_t byte_counts[] = {
        0, 1, 3, 7, 15, 16, 17, 31, 32, 33,
        63, 64, 65, 127, 128, 129, 257, 521
    };
    uint8_t x[524];
    uint8_t y[524];
    uint8_t expected_x[524];
    uint8_t expected_y[524];
    uint8_t output_x[524];
    uint8_t output_y[524];
    uint8_t expected_output_x[524];
    uint8_t expected_output_y[524];

    for (unsigned log = 0; log < 256; ++log)
    {
        for (size_t count_i = 0;
             count_i < sizeof(byte_counts) / sizeof(byte_counts[0]);
             ++count_i)
        {
            const uint64_t bytes = byte_counts[count_i];
            for (size_t i = 0; i < sizeof(x); ++i)
            {
                x[i] = expected_x[i] = static_cast<uint8_t>(
                    i * 73U + log * 11U + count_i);
                y[i] = expected_y[i] = static_cast<uint8_t>(
                    i * 29U + log * 17U + count_i * 3U);
            }
            ReferenceFF8Butterfly2<true>(expected_x + 1, expected_y + 1,
                static_cast<uint8_t>(log), bytes, reference);
            ops.ff8_ifft_butterfly2(
                x + 1, y + 1, static_cast<uint16_t>(log), bytes);
            if (std::memcmp(x, expected_x, sizeof(x)) != 0 ||
                std::memcmp(y, expected_y, sizeof(y)) != 0)
                return false;

            for (size_t i = 0; i < sizeof(x); ++i)
            {
                x[i] = expected_x[i] = static_cast<uint8_t>(
                    i * 61U + log * 7U + count_i);
                y[i] = expected_y[i] = static_cast<uint8_t>(
                    i * 43U + log * 13U + count_i * 5U);
            }
            ReferenceFF8Butterfly2<false>(expected_x + 1, expected_y + 1,
                static_cast<uint8_t>(log), bytes, reference);
            ops.ff8_fft_butterfly2(
                x + 1, y + 1, static_cast<uint16_t>(log), bytes);
            if (std::memcmp(x, expected_x, sizeof(x)) != 0 ||
                std::memcmp(y, expected_y, sizeof(y)) != 0)
                return false;

            for (size_t i = 0; i < sizeof(x); ++i)
            {
                x[i] = static_cast<uint8_t>(
                    i * 37U + log * 19U + count_i);
                y[i] = static_cast<uint8_t>(
                    i * 101U + log * 5U + count_i * 7U);
                output_x[i] = expected_output_x[i] = static_cast<uint8_t>(
                    i * 23U + log + count_i * 11U);
                output_y[i] = expected_output_y[i] = static_cast<uint8_t>(
                    i * 47U + log * 3U + count_i * 13U);
                expected_x[i] = x[i];
                expected_y[i] = y[i];
            }
            ReferenceFF8Butterfly2<true>(expected_x + 1, expected_y + 1,
                static_cast<uint8_t>(log), bytes, reference);
            for (uint64_t i = 0; i < bytes; ++i)
            {
                expected_output_x[i + 1] ^= expected_x[i + 1];
                expected_output_y[i + 1] ^= expected_y[i + 1];
            }
            ops.ff8_ifft_butterfly2_xor(
                x + 1, y + 1, output_x + 1, output_y + 1,
                static_cast<uint16_t>(log), bytes);
            if (std::memcmp(output_x, expected_output_x,
                    sizeof(output_x)) != 0 ||
                std::memcmp(output_y, expected_output_y,
                    sizeof(output_y)) != 0)
                return false;
            // Accumulating butterflies must not consume or modify their input.
            for (size_t i = 0; i < sizeof(x); ++i)
            {
                const uint8_t original_x = static_cast<uint8_t>(
                    i * 37U + log * 19U + count_i);
                const uint8_t original_y = static_cast<uint8_t>(
                    i * 101U + log * 5U + count_i * 7U);
                if (x[i] != original_x || y[i] != original_y)
                    return false;
            }

            for (size_t i = 0; i < sizeof(x); ++i)
            {
                x[i] = expected_x[i] = static_cast<uint8_t>(
                    i * 53U + log * 23U + count_i * 17U);
                y[i] = expected_y[i] = static_cast<uint8_t>(
                    i * 79U + log * 31U + count_i * 19U);
                output_x[i] = expected_output_x[i] = 0xa5;
                output_y[i] = expected_output_y[i] = 0x5a;
            }
            std::memcpy(expected_output_x + 1, x + 1, bytes);
            std::memcpy(expected_output_y + 1, y + 1, bytes);
            if (log == 255)
            {
                for (uint64_t i = 0; i < bytes; ++i)
                    expected_output_y[i + 1] ^= expected_output_x[i + 1];
            }
            else
            {
                ReferenceFF8Butterfly2<false>(
                    expected_output_x + 1, expected_output_y + 1,
                    static_cast<uint8_t>(log), bytes, reference);
            }
            ops.ff8_fft_butterfly2_out(
                x + 1, y + 1, output_x + 1, output_y + 1,
                static_cast<uint16_t>(log), bytes);
            if (std::memcmp(output_x, expected_output_x,
                    sizeof(output_x)) != 0 ||
                std::memcmp(output_y, expected_output_y,
                    sizeof(output_y)) != 0 ||
                std::memcmp(x, expected_x, sizeof(x)) != 0 ||
                std::memcmp(y, expected_y, sizeof(y)) != 0)
                return false;
        }
    }
    return true;
}

template<bool Inverse>
static void ReferenceFF8Butterfly4(
    uint8_t* value0,
    uint8_t* value1,
    uint8_t* value2,
    uint8_t* value3,
    uint16_t log01,
    uint16_t log23,
    uint16_t log02,
    uint64_t byte_count,
    FF8MultiplyLog reference)
{
    static const uint16_t kZeroSkew = 255;
    if (Inverse)
    {
        if (log01 == kZeroSkew)
            for (uint64_t i = 0; i < byte_count; ++i)
                value1[i] ^= value0[i];
        else
            ReferenceFF8Butterfly2<true>(value0, value1,
                static_cast<uint8_t>(log01), byte_count, reference);
        if (log23 == kZeroSkew)
            for (uint64_t i = 0; i < byte_count; ++i)
                value3[i] ^= value2[i];
        else
            ReferenceFF8Butterfly2<true>(value2, value3,
                static_cast<uint8_t>(log23), byte_count, reference);
        if (log02 == kZeroSkew)
            for (uint64_t i = 0; i < byte_count; ++i)
            {
                value2[i] ^= value0[i];
                value3[i] ^= value1[i];
            }
        else
        {
            ReferenceFF8Butterfly2<true>(value0, value2,
                static_cast<uint8_t>(log02), byte_count, reference);
            ReferenceFF8Butterfly2<true>(value1, value3,
                static_cast<uint8_t>(log02), byte_count, reference);
        }
    }
    else
    {
        if (log02 == kZeroSkew)
            for (uint64_t i = 0; i < byte_count; ++i)
            {
                value2[i] ^= value0[i];
                value3[i] ^= value1[i];
            }
        else
        {
            ReferenceFF8Butterfly2<false>(value0, value2,
                static_cast<uint8_t>(log02), byte_count, reference);
            ReferenceFF8Butterfly2<false>(value1, value3,
                static_cast<uint8_t>(log02), byte_count, reference);
        }
        if (log01 == kZeroSkew)
            for (uint64_t i = 0; i < byte_count; ++i)
                value1[i] ^= value0[i];
        else
            ReferenceFF8Butterfly2<false>(value0, value1,
                static_cast<uint8_t>(log01), byte_count, reference);
        if (log23 == kZeroSkew)
            for (uint64_t i = 0; i < byte_count; ++i)
                value3[i] ^= value2[i];
        else
            ReferenceFF8Butterfly2<false>(value2, value3,
                static_cast<uint8_t>(log23), byte_count, reference);
    }
}

static bool TestFF8Butterflies4(const Ops& ops, FF8MultiplyLog reference)
{
    static const uint16_t log_sets[][3] = {
        { 255, 255, 255 },
        { 0, 0, 0 },
        { 1, 2, 3 },
        { 254, 253, 252 },
        { 255, 0, 1 },
        { 2, 255, 3 },
        { 4, 5, 255 },
        { 255, 255, 7 },
        { 255, 11, 255 }
    };
    static const uint64_t byte_counts[] = {
        0, 1, 3, 7, 15, 16, 17, 31, 32, 33,
        63, 64, 65, 127, 128, 129, 257, 521,
        1023, 1024, 1025
    };
    uint8_t values[4][1028];
    uint8_t expected[4][1028];
    uint8_t inputs[4][1028];
    uint8_t original_inputs[4][1028];
    uint8_t outputs[4][1028];
    for (size_t set_i = 0;
         set_i < sizeof(log_sets) / sizeof(log_sets[0]); ++set_i)
        for (size_t count_i = 0;
             count_i < sizeof(byte_counts) / sizeof(byte_counts[0]);
             ++count_i)
        {
            const uint64_t bytes = byte_counts[count_i];
            for (unsigned lane = 0; lane < 4; ++lane)
                for (size_t i = 0; i < sizeof(values[lane]); ++i)
                    values[lane][i] = expected[lane][i] =
                        static_cast<uint8_t>(i * (29U + lane * 12U) +
                            set_i * 31U + count_i * 7U + lane);
            ReferenceFF8Butterfly4<true>(
                expected[0] + 1, expected[1] + 1,
                expected[2] + 1, expected[3] + 1,
                log_sets[set_i][0], log_sets[set_i][1],
                log_sets[set_i][2], bytes, reference);
            ops.ff8_ifft_butterfly4(
                values[0] + 1, values[1] + 1,
                values[2] + 1, values[3] + 1,
                log_sets[set_i][0], log_sets[set_i][1],
                log_sets[set_i][2], bytes);
            if (std::memcmp(values, expected, sizeof(values)) != 0)
                return false;

            for (unsigned lane = 0; lane < 4; ++lane)
                for (size_t i = 0; i < sizeof(values[lane]); ++i)
                    values[lane][i] = expected[lane][i] =
                        static_cast<uint8_t>(i * (43U + lane * 10U) +
                            set_i * 17U + count_i * 11U + lane * 3U);
            ReferenceFF8Butterfly4<false>(
                expected[0] + 1, expected[1] + 1,
                expected[2] + 1, expected[3] + 1,
                log_sets[set_i][0], log_sets[set_i][1],
                log_sets[set_i][2], bytes, reference);
            ops.ff8_fft_butterfly4(
                values[0] + 1, values[1] + 1,
                values[2] + 1, values[3] + 1,
                log_sets[set_i][0], log_sets[set_i][1],
                log_sets[set_i][2], bytes);
            if (std::memcmp(values, expected, sizeof(values)) != 0)
                return false;

            for (unsigned lane = 0; lane < 4; ++lane)
                for (size_t i = 0; i < sizeof(inputs[lane]); ++i)
                {
                    inputs[lane][i] = original_inputs[lane][i] =
                        static_cast<uint8_t>(i * (59U + lane * 14U) +
                            set_i * 37U + count_i * 23U + lane * 5U);
                    outputs[lane][i] = expected[lane][i] =
                        static_cast<uint8_t>(0x91U + lane * 13U);
                }
            for (unsigned lane = 0; lane < 4; ++lane)
                std::memcpy(expected[lane] + 1, inputs[lane] + 1, bytes);
            ReferenceFF8Butterfly4<false>(
                expected[0] + 1, expected[1] + 1,
                expected[2] + 1, expected[3] + 1,
                log_sets[set_i][0], log_sets[set_i][1],
                log_sets[set_i][2], bytes, reference);
            ops.ff8_fft_butterfly4_out(
                inputs[0] + 1, inputs[1] + 1,
                inputs[2] + 1, inputs[3] + 1,
                outputs[0] + 1, outputs[1] + 1,
                outputs[2] + 1, outputs[3] + 1,
                log_sets[set_i][0], log_sets[set_i][1],
                log_sets[set_i][2], bytes);
            if (std::memcmp(outputs, expected, sizeof(outputs)) != 0 ||
                std::memcmp(inputs, original_inputs, sizeof(inputs)) != 0)
                return false;
        }
    return true;
}

static bool TestFF8ButterflyRanges(const Ops& ops)
{
    static const unsigned kDistance = 3;
    static const unsigned kLaneCount = kDistance * 4;
    static const uint16_t log_sets[][3] = {
        { 255, 0, 7 }, { 1, 2, 3 }
    };
    // The leaf KAT already exercises the 1,025-byte cutoff tail.  Keep this
    // range-specific test below common 64-KiB caller-stack budgets while still
    // spanning multiple AVX2 vectors and a byte tail.
    static const uint64_t byte_counts[] = { 17, 129 };
    uint8_t actual[kLaneCount][132];
    uint8_t expected[kLaneCount][132];
    uint8_t xor_actual[kLaneCount][132];
    uint8_t xor_expected[kLaneCount][132];
    void* actual_pointers[kLaneCount];
    void* expected_pointers[kLaneCount];
    void* xor_actual_pointers[kLaneCount];
    for (size_t set_i = 0;
         set_i < sizeof(log_sets) / sizeof(log_sets[0]); ++set_i)
        for (size_t count_i = 0;
             count_i < sizeof(byte_counts) / sizeof(byte_counts[0]);
             ++count_i)
        {
            const uint64_t bytes = byte_counts[count_i];
            for (unsigned lane = 0; lane < kLaneCount; ++lane)
            {
                for (size_t i = 0; i < sizeof(actual[lane]); ++i)
                {
                    actual[lane][i] = expected[lane][i] =
                        static_cast<uint8_t>(
                            i * (17U + lane * 6U) + set_i * 31U + lane);
                    xor_actual[lane][i] = xor_expected[lane][i] =
                        static_cast<uint8_t>(0xa5U + i * 3U + lane * 11U);
                }
                actual_pointers[lane] = actual[lane] + 1;
                expected_pointers[lane] = expected[lane] + 1;
                xor_actual_pointers[lane] = xor_actual[lane] + 1;
            }
            for (unsigned i = 0; i < kDistance; ++i)
                ops.ff8_ifft_butterfly4(
                    expected_pointers[i], expected_pointers[i + kDistance],
                    expected_pointers[i + kDistance * 2U],
                    expected_pointers[i + kDistance * 3U],
                    log_sets[set_i][0], log_sets[set_i][1],
                    log_sets[set_i][2], bytes);
            ops.ff8_ifft_butterfly4_range(
                actual_pointers, kDistance,
                log_sets[set_i][0], log_sets[set_i][1],
                log_sets[set_i][2], bytes, true);
            if (std::memcmp(actual, expected, sizeof(actual)) != 0)
                return false;

            for (unsigned lane = 0; lane < kLaneCount; ++lane)
                for (size_t i = 0; i < sizeof(actual[lane]); ++i)
                    actual[lane][i] = expected[lane][i] =
                        static_cast<uint8_t>(
                            i * (43U + lane * 4U) + set_i * 13U + lane * 3U);
            for (unsigned i = 0; i < kDistance; ++i)
                ops.ff8_fft_butterfly4(
                    expected_pointers[i], expected_pointers[i + kDistance],
                    expected_pointers[i + kDistance * 2U],
                    expected_pointers[i + kDistance * 3U],
                    log_sets[set_i][0], log_sets[set_i][1],
                    log_sets[set_i][2], bytes);
            ops.ff8_fft_butterfly4_range(
                actual_pointers, kDistance,
                log_sets[set_i][0], log_sets[set_i][1],
                log_sets[set_i][2], bytes, true);
            if (std::memcmp(actual, expected, sizeof(actual)) != 0)
                return false;

            for (unsigned lane = 0; lane < kLaneCount; ++lane)
                for (size_t i = 0; i < sizeof(actual[lane]); ++i)
                    actual[lane][i] = expected[lane][i] =
                        static_cast<uint8_t>(
                            i * (61U + lane * 2U) + set_i * 19U + lane * 5U);
            for (unsigned i = 0; i < kDistance; ++i)
                ops.ff8_ifft_butterfly4(
                    expected_pointers[i], expected_pointers[i + kDistance],
                    expected_pointers[i + kDistance * 2U],
                    expected_pointers[i + kDistance * 3U],
                    log_sets[set_i][0], log_sets[set_i][1],
                    log_sets[set_i][2], bytes);
            for (unsigned lane = 0; lane < kLaneCount; ++lane)
                for (uint64_t i = 0; i < bytes; ++i)
                    xor_expected[lane][i + 1] ^= expected[lane][i + 1];
            ops.ff8_ifft_butterfly4_xor_range(
                actual_pointers, xor_actual_pointers, kDistance,
                log_sets[set_i][0], log_sets[set_i][1],
                log_sets[set_i][2], bytes);
            if (std::memcmp(actual, expected, sizeof(actual)) != 0 ||
                std::memcmp(xor_actual, xor_expected,
                    sizeof(xor_actual)) != 0)
                return false;
        }
    return true;
}

#endif // LEO_HAS_FF8

#ifdef LEO_HAS_FF16
static void FillFF16(
    uint8_t* bytes,
    uint64_t byte_count,
    uint32_t seed)
{
    uint64_t offset = 0;
    uint32_t index = seed;
    while (offset <= byte_count && byte_count - offset >= 64)
    {
        for (unsigned i = 0; i < 32; ++i, ++index)
        {
            const uint16_t value = static_cast<uint16_t>(
                index * 40503U + 0x5a17U);
            bytes[offset + i] = static_cast<uint8_t>(value);
            bytes[offset + 32 + i] = static_cast<uint8_t>(value >> 8);
        }
        offset += 64;
    }
    const uint64_t symbols = (byte_count - offset) / 2;
    for (uint64_t i = 0; i < symbols; ++i, ++index)
    {
        const uint16_t value = static_cast<uint16_t>(
            index * 40503U + 0x5a17U);
        bytes[offset + i] = static_cast<uint8_t>(value);
        bytes[offset + symbols + i] = static_cast<uint8_t>(value >> 8);
    }
}

template<bool Add>
static void ReferenceFF16(
    uint8_t* output,
    const uint8_t* input,
    uint16_t log,
    uint64_t byte_count,
    FF16MultiplyLog reference)
{
    uint64_t offset = 0;
    while (byte_count - offset >= 64)
    {
        for (unsigned i = 0; i < 32; ++i)
        {
            const uint16_t value = static_cast<uint16_t>(input[offset + i] |
                (static_cast<unsigned>(input[offset + 32 + i]) << 8));
            const uint16_t product = reference(value, log);
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
    const uint64_t symbols = (byte_count - offset) / 2;
    for (uint64_t i = 0; i < symbols; ++i)
    {
        const uint16_t value = static_cast<uint16_t>(input[offset + i] |
            (static_cast<unsigned>(input[offset + symbols + i]) << 8));
        const uint16_t product = reference(value, log);
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

static bool TestFF16(const Ops& ops, FF16MultiplyLog reference)
{
    static const uint16_t logs[] = {
        0, 1, 2, 3, 15, 16, 255, 256, 4095, 32767, 65534, 65535
    };
    static const uint64_t byte_counts[] = {
        2, 6, 30, 62, 64, 66, 94, 126, 128, 130, 194
    };
    uint8_t source[198];
    uint8_t output[198];
    uint8_t expected[198];

    for (size_t size_i = 0;
         size_i < sizeof(byte_counts) / sizeof(byte_counts[0]);
         ++size_i)
    {
        const uint64_t bytes = byte_counts[size_i];
        std::memset(source, 0x3c, sizeof(source));
        FillFF16(source + 1, bytes, static_cast<uint32_t>(size_i * 41U));
        for (size_t log_i = 0; log_i < sizeof(logs) / sizeof(logs[0]); ++log_i)
        {
            const uint16_t log = logs[log_i];
            std::memset(output, 0xa5, sizeof(output));
            std::memset(expected, 0xa5, sizeof(expected));
            ReferenceFF16<false>(expected + 1, source + 1,
                log, bytes, reference);
            ops.ff16_multiply(output + 1, source + 1, log, bytes);
            if (std::memcmp(output, expected, sizeof(output)) != 0)
                return false;

            for (size_t i = 0; i < sizeof(output); ++i)
                output[i] = expected[i] = static_cast<uint8_t>(
                    (i * 17U + log_i * 13U) & 255U);
            ReferenceFF16<true>(expected + 1, source + 1,
                log, bytes, reference);
            ops.ff16_multiply_add(output + 1, source + 1, log, bytes);
            if (std::memcmp(output, expected, sizeof(output)) != 0)
                return false;
        }
    }
    return true;
}

template<bool Inverse>
static void ReferenceFF16Butterfly2(
    uint8_t* x,
    uint8_t* y,
    uint16_t log,
    uint64_t byte_count,
    FF16MultiplyLog reference)
{
    uint8_t product[194];
    if (Inverse)
        for (uint64_t i = 0; i < byte_count; ++i)
            y[i] ^= x[i];
    ReferenceFF16<false>(product, y, log, byte_count, reference);
    for (uint64_t i = 0; i < byte_count; ++i)
        x[i] ^= product[i];
    if (!Inverse)
        for (uint64_t i = 0; i < byte_count; ++i)
            y[i] ^= x[i];
}

static bool TestFF16Butterflies(const Ops& ops, FF16MultiplyLog reference)
{
    static const uint16_t logs[] = {
        0, 1, 2, 3, 15, 16, 255, 256, 4095, 32767, 65534, 65535
    };
    static const uint64_t byte_counts[] = {
        0, 2, 6, 30, 62, 64, 66, 94, 126, 128, 130, 194
    };
    uint8_t x[198];
    uint8_t y[198];
    uint8_t expected_x[198];
    uint8_t expected_y[198];
    uint8_t output_x[198];
    uint8_t output_y[198];
    uint8_t expected_output_x[198];
    uint8_t expected_output_y[198];
    for (size_t size_i = 0;
         size_i < sizeof(byte_counts) / sizeof(byte_counts[0]);
         ++size_i)
    {
        const uint64_t bytes = byte_counts[size_i];
        for (size_t log_i = 0; log_i < sizeof(logs) / sizeof(logs[0]); ++log_i)
        {
            const uint16_t log = logs[log_i];
            std::memset(x, 0xa5, sizeof(x));
            std::memset(y, 0x5a, sizeof(y));
            FillFF16(x + 1, bytes,
                static_cast<uint32_t>(size_i * 41U + log_i));
            FillFF16(y + 1, bytes,
                static_cast<uint32_t>(size_i * 67U + log_i * 3U));
            std::memcpy(expected_x, x, sizeof(x));
            std::memcpy(expected_y, y, sizeof(y));
            ReferenceFF16Butterfly2<true>(expected_x + 1, expected_y + 1,
                log, bytes, reference);
            ops.ff16_ifft_butterfly2(x + 1, y + 1, log, bytes);
            if (std::memcmp(x, expected_x, sizeof(x)) != 0 ||
                std::memcmp(y, expected_y, sizeof(y)) != 0)
                return false;

            std::memset(x, 0x3c, sizeof(x));
            std::memset(y, 0xc3, sizeof(y));
            FillFF16(x + 1, bytes,
                static_cast<uint32_t>(size_i * 73U + log_i * 5U));
            FillFF16(y + 1, bytes,
                static_cast<uint32_t>(size_i * 89U + log_i * 7U));
            std::memcpy(expected_x, x, sizeof(x));
            std::memcpy(expected_y, y, sizeof(y));
            ReferenceFF16Butterfly2<false>(expected_x + 1, expected_y + 1,
                log, bytes, reference);
            ops.ff16_fft_butterfly2(x + 1, y + 1, log, bytes);
            if (std::memcmp(x, expected_x, sizeof(x)) != 0 ||
                std::memcmp(y, expected_y, sizeof(y)) != 0)
                return false;

            std::memset(x, 0x96, sizeof(x));
            std::memset(y, 0x69, sizeof(y));
            FillFF16(x + 1, bytes,
                static_cast<uint32_t>(size_i * 97U + log_i * 11U));
            FillFF16(y + 1, bytes,
                static_cast<uint32_t>(size_i * 103U + log_i * 13U));
            std::memcpy(expected_x, x, sizeof(x));
            std::memcpy(expected_y, y, sizeof(y));
            std::memset(output_x, 0xa5, sizeof(output_x));
            std::memset(output_y, 0x5a, sizeof(output_y));
            std::memcpy(expected_output_x, output_x, sizeof(output_x));
            std::memcpy(expected_output_y, output_y, sizeof(output_y));
            std::memcpy(expected_output_x + 1, x + 1, bytes);
            std::memcpy(expected_output_y + 1, y + 1, bytes);
            if (log == 65535)
            {
                for (uint64_t i = 0; i < bytes; ++i)
                    expected_output_y[i + 1] ^= expected_output_x[i + 1];
            }
            else
            {
                ReferenceFF16Butterfly2<false>(
                    expected_output_x + 1, expected_output_y + 1,
                    log, bytes, reference);
            }
            ops.ff16_fft_butterfly2_out(
                x + 1, y + 1, output_x + 1, output_y + 1, log, bytes);
            if (std::memcmp(output_x, expected_output_x,
                    sizeof(output_x)) != 0 ||
                std::memcmp(output_y, expected_output_y,
                    sizeof(output_y)) != 0 ||
                std::memcmp(x, expected_x, sizeof(x)) != 0 ||
                std::memcmp(y, expected_y, sizeof(y)) != 0)
                return false;
        }
    }
    return true;
}

template<bool Inverse>
static void ReferenceFF16Butterfly4(
    uint8_t* value0,
    uint8_t* value1,
    uint8_t* value2,
    uint8_t* value3,
    uint16_t log01,
    uint16_t log23,
    uint16_t log02,
    uint64_t byte_count,
    FF16MultiplyLog reference)
{
    static const uint16_t kZeroSkew = 65535;
    if (Inverse)
    {
        if (log01 == kZeroSkew)
            for (uint64_t i = 0; i < byte_count; ++i)
                value1[i] ^= value0[i];
        else
            ReferenceFF16Butterfly2<true>(
                value0, value1, log01, byte_count, reference);
        if (log23 == kZeroSkew)
            for (uint64_t i = 0; i < byte_count; ++i)
                value3[i] ^= value2[i];
        else
            ReferenceFF16Butterfly2<true>(
                value2, value3, log23, byte_count, reference);
        if (log02 == kZeroSkew)
            for (uint64_t i = 0; i < byte_count; ++i)
            {
                value2[i] ^= value0[i];
                value3[i] ^= value1[i];
            }
        else
        {
            ReferenceFF16Butterfly2<true>(
                value0, value2, log02, byte_count, reference);
            ReferenceFF16Butterfly2<true>(
                value1, value3, log02, byte_count, reference);
        }
    }
    else
    {
        if (log02 == kZeroSkew)
            for (uint64_t i = 0; i < byte_count; ++i)
            {
                value2[i] ^= value0[i];
                value3[i] ^= value1[i];
            }
        else
        {
            ReferenceFF16Butterfly2<false>(
                value0, value2, log02, byte_count, reference);
            ReferenceFF16Butterfly2<false>(
                value1, value3, log02, byte_count, reference);
        }
        if (log01 == kZeroSkew)
            for (uint64_t i = 0; i < byte_count; ++i)
                value1[i] ^= value0[i];
        else
            ReferenceFF16Butterfly2<false>(
                value0, value1, log01, byte_count, reference);
        if (log23 == kZeroSkew)
            for (uint64_t i = 0; i < byte_count; ++i)
                value3[i] ^= value2[i];
        else
            ReferenceFF16Butterfly2<false>(
                value2, value3, log23, byte_count, reference);
    }
}

static bool TestFF16Butterflies4(const Ops& ops, FF16MultiplyLog reference)
{
    static const uint16_t log_sets[][3] = {
        { 65535, 65535, 65535 },
        { 0, 0, 0 },
        { 1, 256, 4095 },
        { 65534, 32767, 255 },
        { 65535, 0, 1 },
        { 2, 65535, 3 },
        { 4, 5, 65535 },
        { 65535, 65535, 32767 },
        { 65535, 16, 65535 }
    };
    static const uint64_t byte_counts[] = {
        0, 2, 6, 30, 62, 64, 66, 94, 126, 128, 130, 194
    };
    uint8_t values[4][198];
    uint8_t expected[4][198];
    uint8_t inputs[4][198];
    uint8_t original_inputs[4][198];
    uint8_t outputs[4][198];
    for (size_t set_i = 0;
         set_i < sizeof(log_sets) / sizeof(log_sets[0]); ++set_i)
        for (size_t count_i = 0;
             count_i < sizeof(byte_counts) / sizeof(byte_counts[0]);
             ++count_i)
        {
            const uint64_t bytes = byte_counts[count_i];
            for (unsigned lane = 0; lane < 4; ++lane)
            {
                std::memset(values[lane], 0xa5U + lane,
                    sizeof(values[lane]));
                FillFF16(values[lane] + 1, bytes,
                    static_cast<uint32_t>(
                        set_i * 101U + count_i * 17U + lane * 7U));
                std::memcpy(expected[lane], values[lane],
                    sizeof(values[lane]));
            }
            ReferenceFF16Butterfly4<true>(
                expected[0] + 1, expected[1] + 1,
                expected[2] + 1, expected[3] + 1,
                log_sets[set_i][0], log_sets[set_i][1],
                log_sets[set_i][2], bytes, reference);
            ops.ff16_ifft_butterfly4(
                values[0] + 1, values[1] + 1,
                values[2] + 1, values[3] + 1,
                log_sets[set_i][0], log_sets[set_i][1],
                log_sets[set_i][2], bytes);
            if (std::memcmp(values, expected, sizeof(values)) != 0)
                return false;

            for (unsigned lane = 0; lane < 4; ++lane)
            {
                std::memset(values[lane], 0x3cU + lane,
                    sizeof(values[lane]));
                FillFF16(values[lane] + 1, bytes,
                    static_cast<uint32_t>(
                        set_i * 149U + count_i * 29U + lane * 11U));
                std::memcpy(expected[lane], values[lane],
                    sizeof(values[lane]));
            }
            ReferenceFF16Butterfly4<false>(
                expected[0] + 1, expected[1] + 1,
                expected[2] + 1, expected[3] + 1,
                log_sets[set_i][0], log_sets[set_i][1],
                log_sets[set_i][2], bytes, reference);
            ops.ff16_fft_butterfly4(
                values[0] + 1, values[1] + 1,
                values[2] + 1, values[3] + 1,
                log_sets[set_i][0], log_sets[set_i][1],
                log_sets[set_i][2], bytes);
            if (std::memcmp(values, expected, sizeof(values)) != 0)
                return false;

            for (unsigned lane = 0; lane < 4; ++lane)
            {
                std::memset(inputs[lane], 0x71U + lane, sizeof(inputs[lane]));
                FillFF16(inputs[lane] + 1, bytes,
                    static_cast<uint32_t>(
                        set_i * 173U + count_i * 31U + lane * 19U));
                std::memcpy(original_inputs[lane], inputs[lane],
                    sizeof(inputs[lane]));
                std::memset(outputs[lane], 0x91U + lane,
                    sizeof(outputs[lane]));
                std::memcpy(expected[lane], outputs[lane],
                    sizeof(outputs[lane]));
                std::memcpy(expected[lane] + 1, inputs[lane] + 1, bytes);
            }
            ReferenceFF16Butterfly4<false>(
                expected[0] + 1, expected[1] + 1,
                expected[2] + 1, expected[3] + 1,
                log_sets[set_i][0], log_sets[set_i][1],
                log_sets[set_i][2], bytes, reference);
            ops.ff16_fft_butterfly4_out(
                inputs[0] + 1, inputs[1] + 1,
                inputs[2] + 1, inputs[3] + 1,
                outputs[0] + 1, outputs[1] + 1,
                outputs[2] + 1, outputs[3] + 1,
                log_sets[set_i][0], log_sets[set_i][1],
                log_sets[set_i][2], bytes);
            if (std::memcmp(outputs, expected, sizeof(outputs)) != 0 ||
                std::memcmp(inputs, original_inputs, sizeof(inputs)) != 0)
                return false;
        }
    return true;
}

static bool TestFF16ButterflyRanges(const Ops& ops)
{
    static const unsigned kDistance = 3;
    static const unsigned kLaneCount = kDistance * 4;
    static const uint16_t log_sets[][3] = {
        { 65535, 0, 17 }, { 1, 256, 4095 }
    };
    static const uint64_t byte_counts[] = { 64, 128, 130 };
    uint8_t actual[kLaneCount][198];
    uint8_t expected[kLaneCount][198];
    void* actual_pointers[kLaneCount];
    void* expected_pointers[kLaneCount];
    for (size_t set_i = 0;
         set_i < sizeof(log_sets) / sizeof(log_sets[0]); ++set_i)
        for (size_t count_i = 0;
             count_i < sizeof(byte_counts) / sizeof(byte_counts[0]);
             ++count_i)
            for (unsigned inverse = 0; inverse < 2; ++inverse)
                for (unsigned fused = 0; fused < 2; ++fused)
                {
                    const uint64_t bytes = byte_counts[count_i];
                    for (unsigned lane = 0; lane < kLaneCount; ++lane)
                    {
                        for (size_t i = 0; i < sizeof(actual[lane]); ++i)
                            actual[lane][i] = expected[lane][i] =
                                static_cast<uint8_t>(
                                    i * (23U + lane * 4U) + set_i * 29U +
                                    inverse * 37U + fused * 41U + lane);
                        actual_pointers[lane] = actual[lane] + 1;
                        expected_pointers[lane] = expected[lane] + 1;
                    }
                    for (unsigned i = 0; i < kDistance; ++i)
                    {
                        if (inverse)
                            ops.ff16_ifft_butterfly4(
                                expected_pointers[i],
                                expected_pointers[i + kDistance],
                                expected_pointers[i + kDistance * 2U],
                                expected_pointers[i + kDistance * 3U],
                                log_sets[set_i][0], log_sets[set_i][1],
                                log_sets[set_i][2], bytes);
                        else
                            ops.ff16_fft_butterfly4(
                                expected_pointers[i],
                                expected_pointers[i + kDistance],
                                expected_pointers[i + kDistance * 2U],
                                expected_pointers[i + kDistance * 3U],
                                log_sets[set_i][0], log_sets[set_i][1],
                                log_sets[set_i][2], bytes);
                    }
                    (inverse ? ops.ff16_ifft_butterfly4_range :
                               ops.ff16_fft_butterfly4_range)(
                        actual_pointers, kDistance,
                        log_sets[set_i][0], log_sets[set_i][1],
                        log_sets[set_i][2], bytes, fused != 0);
                    if (std::memcmp(actual, expected, sizeof(actual)) != 0)
                        return false;
                }
    return true;
}

#endif // LEO_HAS_FF16

static bool TestXor(const Ops& ops)
{
    static const uint64_t byte_counts[] = {
        0, 1, 3, 7, 15, 16, 17, 31, 32, 33,
        63, 64, 65, 127, 128, 129, 257
    };
    uint8_t source[260];
    uint8_t output[260];
    uint8_t expected[260];
    uint8_t sources4[3][260];
    uint8_t outputs4[4][260];
    uint8_t expected4[4][260];
    for (size_t i = 0; i < sizeof(source); ++i)
    {
        source[i] = static_cast<uint8_t>((i * 101U + 7U) & 255U);
        for (unsigned source_i = 0; source_i < 3; ++source_i)
            sources4[source_i][i] = static_cast<uint8_t>(
                i * (37U + source_i * 18U) + source_i * 53U + 11U);
    }
    for (size_t count_i = 0;
         count_i < sizeof(byte_counts) / sizeof(byte_counts[0]);
         ++count_i)
    {
        for (size_t i = 0; i < sizeof(output); ++i)
            output[i] = expected[i] =
                static_cast<uint8_t>((i * 31U + count_i) & 255U);
        const uint64_t bytes = byte_counts[count_i];
        for (uint64_t i = 0; i < bytes; ++i)
            expected[i + 1] ^= source[i + 1];
        ops.xor_memory(output + 1, source + 1, bytes);
        if (std::memcmp(output, expected, sizeof(output)) != 0)
            return false;

        for (size_t i = 0; i < sizeof(output); ++i)
            output[i] = expected[i] =
                static_cast<uint8_t>((i * 43U + count_i * 3U) & 255U);
        for (uint64_t i = 0; i < bytes; ++i)
            expected[i + 1] ^= source[i + 1] ^ sources4[0][i + 1];
        ops.xor_memory_2to1(
            output + 1, source + 1, sources4[0] + 1, bytes);
        if (std::memcmp(output, expected, sizeof(output)) != 0)
            return false;

        // The two read-only inputs may be the same range.  Their contribution
        // cancels exactly and the destination, including guards, is unchanged.
        std::memcpy(expected, output, sizeof(output));
        ops.xor_memory_2to1(
            output + 1, source + 1, source + 1, bytes);
        if (std::memcmp(output, expected, sizeof(output)) != 0)
            return false;

        for (unsigned pair = 0; pair < 4; ++pair)
            for (size_t i = 0; i < sizeof(outputs4[pair]); ++i)
                outputs4[pair][i] = expected4[pair][i] =
                    static_cast<uint8_t>(
                        i * (29U + pair * 6U) + count_i * 7U + pair);
        // Pair 2 deliberately shares pair 0's read-only input.  Public shard
        // inputs may alias each other even though every destination is
        // disjoint from every source and other destination.
        const unsigned source_indices[4] = { 0, 1, 0, 2 };
        for (unsigned pair = 0; pair < 4; ++pair)
            for (uint64_t i = 0; i < bytes; ++i)
                expected4[pair][i + 1] ^=
                    sources4[source_indices[pair]][i + 1];
        ops.xor_memory4(
            outputs4[0] + 1, sources4[0] + 1,
            outputs4[1] + 1, sources4[1] + 1,
            outputs4[2] + 1, sources4[0] + 1,
            outputs4[3] + 1, sources4[2] + 1, bytes);
        if (std::memcmp(outputs4, expected4, sizeof(outputs4)) != 0)
            return false;
    }
    return true;
}

static bool TestOps(const Ops& ops, const InitializeArgs& args)
{
    if (!ops.name || !ops.xor_memory || !ops.xor_memory_2to1 ||
        !ops.xor_memory4 || !TestXor(ops))
        return false;
#ifdef LEO_HAS_FF8
    if (!args.ff8_multiply_log || !ops.ff8_multiply ||
        !ops.ff8_multiply_add || !ops.ff8_ifft_butterfly2 ||
        !ops.ff8_fft_butterfly2 || !ops.ff8_fft_butterfly2_out ||
        !ops.ff8_ifft_butterfly2_xor || !ops.ff8_ifft_butterfly4 ||
        !ops.ff8_fft_butterfly4 || !ops.ff8_fft_butterfly4_out ||
        !ops.ff8_ifft_butterfly4_range ||
        !ops.ff8_fft_butterfly4_range ||
        !ops.ff8_ifft_butterfly4_xor_range ||
        !TestFF8(ops, args.ff8_multiply_log) ||
        !TestFF8Butterflies(ops, args.ff8_multiply_log) ||
        !TestFF8Butterflies4(ops, args.ff8_multiply_log) ||
        !TestFF8ButterflyRanges(ops))
        return false;
#else
    if (ops.ff8_multiply || ops.ff8_multiply_add ||
        ops.ff8_ifft_butterfly2 || ops.ff8_fft_butterfly2 ||
        ops.ff8_fft_butterfly2_out || ops.ff8_ifft_butterfly2_xor ||
        ops.ff8_ifft_butterfly4 || ops.ff8_fft_butterfly4 ||
        ops.ff8_fft_butterfly4_out || ops.ff8_ifft_butterfly4_range ||
        ops.ff8_fft_butterfly4_range ||
        ops.ff8_ifft_butterfly4_xor_range)
        return false;
#endif
#ifdef LEO_HAS_FF16
    if (!args.ff16_multiply_log || !ops.ff16_multiply ||
        !ops.ff16_multiply_add || !ops.ff16_ifft_butterfly2 ||
        !ops.ff16_fft_butterfly2 || !ops.ff16_fft_butterfly2_out ||
        !ops.ff16_ifft_butterfly4 || !ops.ff16_fft_butterfly4 ||
        !ops.ff16_fft_butterfly4_out ||
        !ops.ff16_ifft_butterfly4_range ||
        !ops.ff16_fft_butterfly4_range ||
        !TestFF16(ops, args.ff16_multiply_log) ||
        !TestFF16Butterflies(ops, args.ff16_multiply_log) ||
        !TestFF16Butterflies4(ops, args.ff16_multiply_log) ||
        !TestFF16ButterflyRanges(ops))
        return false;
#else
    if (ops.ff16_multiply || ops.ff16_multiply_add ||
        ops.ff16_ifft_butterfly2 || ops.ff16_fft_butterfly2 ||
        ops.ff16_fft_butterfly2_out || ops.ff16_ifft_butterfly4 ||
        ops.ff16_fft_butterfly4 || ops.ff16_fft_butterfly4_out ||
        ops.ff16_ifft_butterfly4_range ||
        ops.ff16_fft_butterfly4_range)
        return false;
#endif
    return true;
}

static std::mutex& GetQualificationMutex()
{
    static std::mutex mutex;
    return mutex;
}

static bool RegisterQualifiedOps(
    leo2_backend expected_kind,
    const Ops* ops,
    const InitializeArgs& args)
{
    if (!ops || ops->kind != expected_kind)
        return false;
#ifdef LEO2_ENABLE_TEST_HOOKS
    if (ConsumeTestFault(KATFaultFor(expected_kind)))
        return false;
#endif
    if (!TestOps(*ops, args))
        return false;
    QualifiedOps[expected_kind] = ops;
    return true;
}

static leo2_backend SelectBackend(const X86Features& features)
{
    // Some diagnostic-only builds compile this dispatcher without any
    // separately linked x86 backend.  Keep those reviewed configurations
    // warning-clean while the feature-dependent branches compile away.
    (void)features;
    leo2_backend selected_kind = LEO2_BACKEND_SCALAR;

#if defined(LEO2_BACKEND_FORCE_SCALAR)
    selected_kind = LEO2_BACKEND_SCALAR;
#elif defined(LEO2_BACKEND_FORCE_SSSE3)
# if defined(LEO2_HAVE_SSSE3_BACKEND)
    if (!features.ssse3)
        return LEO2_BACKEND_AUTO;
    selected_kind = LEO2_BACKEND_SSSE3;
# else
    return LEO2_BACKEND_AUTO;
# endif
#elif defined(LEO2_BACKEND_FORCE_AVX2)
# if defined(LEO2_HAVE_AVX2_BACKEND)
    if (!features.avx2)
        return LEO2_BACKEND_AUTO;
    selected_kind = LEO2_BACKEND_AVX2;
# else
    return LEO2_BACKEND_AUTO;
# endif
#else
# if defined(LEO2_HAVE_AVX2_BACKEND)
    if (features.avx2)
        selected_kind = LEO2_BACKEND_AVX2;
# endif
# if defined(LEO2_HAVE_SSSE3_BACKEND)
    if (selected_kind == LEO2_BACKEND_SCALAR && features.ssse3)
        selected_kind = LEO2_BACKEND_SSSE3;
# endif
#endif
    return selected_kind;
}

bool Initialize(const InitializeArgs& args)
{
    if (SelectedOps)
        return SelfTestPassed;

    const X86Features features = DetectX86Features();
    (void)features;
    const leo2_backend selected_kind = SelectBackend(features);
    if (selected_kind == LEO2_BACKEND_AUTO)
        return false;

    // The selected process default is mandatory.  Lower tables are initialized
    // lazily by an explicit context request, so legacy/AUTO-only applications
    // do not pay their allocation and table-generation costs.
    const Ops* selected_ops = NULL;
    if (selected_kind == LEO2_BACKEND_SCALAR)
        selected_ops = InitializeScalar(args);

#if defined(LEO2_HAVE_SSSE3_BACKEND)
    if (selected_kind == LEO2_BACKEND_SSSE3)
        selected_ops = InitializeSSSE3(args);
#endif

#if defined(LEO2_HAVE_AVX2_BACKEND)
    if (selected_kind == LEO2_BACKEND_AVX2)
        selected_ops = InitializeAVX2(args);
#endif

    if (!selected_ops)
    {
        QualificationStates[selected_kind] = QualificationFailed;
        QualificationFailures[selected_kind] = QualificationOutOfMemory;
        StartupFailure = QualificationOutOfMemory;
        return false;
    }
    if (!RegisterQualifiedOps(selected_kind, selected_ops, args))
    {
        QualificationStates[selected_kind] = QualificationFailed;
        QualificationFailures[selected_kind] = QualificationSelfTestFailed;
        StartupFailure = QualificationSelfTestFailed;
        return false;
    }

    uint32_t qualifiable_mask = 1U << LEO2_BACKEND_SCALAR;
#if defined(LEO2_HAVE_SSSE3_BACKEND)
    if (features.ssse3 && selected_kind >= LEO2_BACKEND_SSSE3)
        qualifiable_mask |= 1U << LEO2_BACKEND_SSSE3;
#endif
#if defined(LEO2_HAVE_AVX2_BACKEND)
    if (features.avx2 && selected_kind >= LEO2_BACKEND_AVX2)
        qualifiable_mask |= 1U << LEO2_BACKEND_AVX2;
#endif

    {
        std::lock_guard<std::mutex> lock(GetQualificationMutex());
        SavedInitializeArgs = args;
        QualifiableBackendMask = qualifiable_mask;
        QualificationStates[selected_kind] = QualificationPassed;
        StartupFailure = QualificationAvailable;
        SelfTestPassed = true;
        SelectedOps = QualifiedOps[selected_kind];
    }
    return true;
}

const Ops& GetOps()
{
    return *SelectedOps;
}

const Ops& GetDefaultOps()
{
    return *SelectedOps;
}

const Ops* GetQualifiedOps(
    leo2_backend requested,
    QualificationStatus* status)
{
    if (status)
        *status = QualificationAvailable;
    if (requested == LEO2_BACKEND_AUTO)
        return SelectedOps;
    const unsigned index = static_cast<unsigned>(requested);
    if (index > static_cast<unsigned>(LEO2_BACKEND_NEON))
    {
        if (status)
            *status = QualificationUnavailable;
        return NULL;
    }

    std::lock_guard<std::mutex> lock(GetQualificationMutex());
    if ((QualifiableBackendMask & (1U << index)) == 0)
    {
        if (status)
            *status = QualificationUnavailable;
        return NULL;
    }
    if (QualificationStates[index] == QualificationPassed)
        return QualifiedOps[index];
    if (QualificationStates[index] == QualificationFailed)
    {
        if (status)
            *status = QualificationFailures[index];
        return NULL;
    }

    const Ops* candidate = NULL;
    switch (requested)
    {
    case LEO2_BACKEND_SCALAR:
        candidate = InitializeScalar(SavedInitializeArgs);
        break;
#if defined(LEO2_HAVE_SSSE3_BACKEND)
    case LEO2_BACKEND_SSSE3:
        candidate = InitializeSSSE3(SavedInitializeArgs);
        break;
#endif
#if defined(LEO2_HAVE_AVX2_BACKEND)
    case LEO2_BACKEND_AVX2:
        candidate = InitializeAVX2(SavedInitializeArgs);
        break;
#endif
    default:
        break;
    }
    if (!candidate)
    {
        QualificationStates[index] = QualificationFailed;
        QualificationFailures[index] = QualificationOutOfMemory;
        if (status)
            *status = QualificationOutOfMemory;
        return NULL;
    }
    if (!RegisterQualifiedOps(requested, candidate, SavedInitializeArgs))
    {
        QualificationStates[index] = QualificationFailed;
        QualificationFailures[index] = QualificationSelfTestFailed;
        if (status)
            *status = QualificationSelfTestFailed;
        return NULL;
    }
    QualificationStates[index] = QualificationPassed;
    return QualifiedOps[index];
}

leo2_backend SelectedBackend()
{
    return SelectedOps ? SelectedOps->kind : LEO2_BACKEND_AUTO;
}

bool StartupSelfTestPassed()
{
    return SelfTestPassed;
}

QualificationStatus StartupQualificationFailure()
{
    return StartupFailure;
}

#ifdef LEO2_ENABLE_TEST_HOOKS
void TestSetSetupFault(TestSetupFault fault)
{
    TestFaultConsumptions.store(0, std::memory_order_relaxed);
    TestFault.store(static_cast<unsigned>(fault), std::memory_order_release);
}

bool TestSetupFaultPending()
{
    return TestFault.load(std::memory_order_acquire) !=
        static_cast<unsigned>(TestSetupFaultNone);
}

unsigned TestSetupFaultConsumptions()
{
    return TestFaultConsumptions.load(std::memory_order_relaxed);
}

bool TestShouldFailAllocation(leo2_backend backend, bool ff16)
{
    return ConsumeTestFault(AllocationFaultFor(backend, ff16));
}

leo2_backend TestDefaultBackendForHost()
{
    return SelectBackend(DetectX86Features());
}

bool TestGetBackendState(leo2_backend backend, TestBackendState* state)
{
    if (!state)
        return false;
    *state = TestBackendState();
    switch (backend)
    {
    case LEO2_BACKEND_SCALAR:
        TestGetScalarTableState(state);
        break;
# if defined(LEO2_HAVE_SSSE3_BACKEND)
    case LEO2_BACKEND_SSSE3:
        TestGetSSSE3TableState(state);
        break;
# endif
# if defined(LEO2_HAVE_AVX2_BACKEND)
    case LEO2_BACKEND_AVX2:
        TestGetAVX2TableState(state);
        break;
# endif
    default:
        return false;
    }
    const unsigned index = static_cast<unsigned>(backend);
    std::lock_guard<std::mutex> lock(GetQualificationMutex());
    state->qualified = QualifiedOps[index] != NULL &&
        QualificationStates[index] == QualificationPassed;
    if (QualificationStates[index] == QualificationFailed)
        state->failure = QualificationFailures[index];
    else
        state->failure = QualificationAvailable;
    return true;
}
#endif

}} // namespace leopard::backend
