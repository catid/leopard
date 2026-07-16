/*
    Copyright (c) 2017 Christopher A. Taylor.  All rights reserved.
    See Leopard2Backend.h for the full BSD license notice.
*/

#include "Leopard2Backend.h"

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
        }
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
        ops.ff8_ifft_butterfly2 && ops.ff8_fft_butterfly2 &&
        ops.ff8_ifft_butterfly2_xor && ops.ff16_ifft_butterfly2 &&
        ops.ff16_fft_butterfly2 &&
        TestFF8(ops, args.ff8_multiply_log) &&
        TestFF8Butterflies(ops, args.ff8_multiply_log) &&
        TestFF16(ops, args.ff16_multiply_log) &&
        TestFF16Butterflies(ops, args.ff16_multiply_log) && TestXor(ops);
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
    if (!ops || ops->kind != expected_kind || !TestOps(*ops, args))
        return false;
    QualifiedOps[expected_kind] = ops;
    return true;
}

bool Initialize(const InitializeArgs& args)
{
    if (SelectedOps)
        return SelfTestPassed;

    const X86Features features = DetectX86Features();
    (void)features;
    leo2_backend selected_kind = LEO2_BACKEND_SCALAR;

#if defined(LEO2_BACKEND_FORCE_SCALAR)
    selected_kind = LEO2_BACKEND_SCALAR;
#elif defined(LEO2_BACKEND_FORCE_SSSE3)
# if defined(LEO2_HAVE_SSSE3_BACKEND)
    if (!features.ssse3)
        return false;
    selected_kind = LEO2_BACKEND_SSSE3;
# else
    return false;
# endif
#elif defined(LEO2_BACKEND_FORCE_AVX2)
# if defined(LEO2_HAVE_AVX2_BACKEND)
    if (!features.avx2)
        return false;
    selected_kind = LEO2_BACKEND_AVX2;
# else
    return false;
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

    if (!RegisterQualifiedOps(selected_kind, selected_ops, args))
        return false;

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

}} // namespace leopard::backend
