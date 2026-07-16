/*
    Copyright (c) 2017 Christopher A. Taylor.  All rights reserved.
    See Leopard2Backend.h for the full BSD license notice.
*/

#include "Leopard2Backend.h"

#include <cstring>

namespace leopard { namespace backend {

static const Ops* SelectedOps = NULL;
static bool SelfTestPassed = false;

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

static void FillFF16(
    uint8_t* bytes,
    uint64_t byte_count,
    uint32_t seed)
{
    uint64_t offset = 0;
    uint32_t index = seed;
    while (byte_count - offset >= 64)
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

static bool TestXor(const Ops& ops)
{
    static const uint64_t byte_counts[] = {
        0, 1, 3, 7, 15, 16, 17, 31, 32, 33,
        63, 64, 65, 127, 128, 129, 257
    };
    uint8_t source[260];
    uint8_t output[260];
    uint8_t expected[260];
    for (size_t i = 0; i < sizeof(source); ++i)
        source[i] = static_cast<uint8_t>((i * 101U + 7U) & 255U);
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
    }
    return true;
}

static bool TestOps(const Ops& ops, const InitializeArgs& args)
{
    return ops.name && ops.ff8_multiply && ops.ff8_multiply_add &&
        ops.ff16_multiply && ops.ff16_multiply_add && ops.xor_memory &&
        TestFF8(ops, args.ff8_multiply_log) &&
        TestFF16(ops, args.ff16_multiply_log) && TestXor(ops);
}

bool Initialize(const InitializeArgs& args)
{
    if (SelectedOps)
        return SelfTestPassed;

    const X86Features features = DetectX86Features();
    (void)features;
    const Ops* candidate = NULL;

#if defined(LEO2_BACKEND_FORCE_SCALAR)
    candidate = InitializeScalar(args);
#elif defined(LEO2_BACKEND_FORCE_SSSE3)
# if defined(LEO2_HAVE_SSSE3_BACKEND)
    if (features.ssse3)
        candidate = InitializeSSSE3(args);
# endif
#elif defined(LEO2_BACKEND_FORCE_AVX2)
# if defined(LEO2_HAVE_AVX2_BACKEND)
    if (features.avx2)
        candidate = InitializeAVX2(args);
# endif
#else
# if defined(LEO2_HAVE_AVX2_BACKEND)
    if (features.avx2)
        candidate = InitializeAVX2(args);
# endif
# if defined(LEO2_HAVE_SSSE3_BACKEND)
    if (!candidate && features.ssse3)
        candidate = InitializeSSSE3(args);
# endif
    if (!candidate)
        candidate = InitializeScalar(args);
#endif

    if (!candidate || !TestOps(*candidate, args))
        return false;
    SelfTestPassed = true;
    SelectedOps = candidate;
    return true;
}

const Ops& GetOps()
{
    return *SelectedOps;
}

leo2_backend SelectedBackend()
{
    return SelectedOps ? SelectedOps->kind : LEO2_BACKEND_AUTO;
}

bool StartupSelfTestPassed()
{
    return SelfTestPassed;
}

}} // namespace leopard::backend
