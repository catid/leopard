#include "../LeopardFF8Xor.h"
#include "../LeopardFF8XorAVX512Four.h"
#include "../leopard_ff8xor.h"

#include <stdint.h>
#include <stdio.h>
#include <string.h>

#include <array>
#include <vector>

namespace {

typedef std::vector<uint8_t> Buffer;
typedef std::array<Buffer, 4> FourBuffers;

static uint32_t RandomState = UINT32_C(0x4f8c17a5);

static uint8_t NextByte()
{
    RandomState ^= RandomState << 13;
    RandomState ^= RandomState >> 17;
    RandomState ^= RandomState << 5;
    return static_cast<uint8_t>(RandomState >> 24);
}

static void Fill(Buffer& buffer)
{
    for (size_t i = 0; i < buffer.size(); ++i)
        buffer[i] = NextByte();
}

static bool EqualFour(
    const FourBuffers& expected,
    const FourBuffers& actual,
    unsigned tuple_index,
    bool inverse,
    leopard::ff8xor::FourBufferMode mode)
{
    for (unsigned index = 0; index < 4; ++index)
    {
        if (expected[index] != actual[index])
        {
            fprintf(stderr,
                "four-buffer mismatch: tuple=%u inverse=%u mode=%u "
                "buffer=%u\n",
                tuple_index,
                inverse ? 1U : 0U,
                static_cast<unsigned>(mode),
                index);
            return false;
        }
    }
    return true;
}

static void ApplyReference(
    FourBuffers& buffers,
    bool inverse,
    uint8_t skew01,
    uint8_t skew23,
    uint8_t skew02)
{
    const uint64_t bytes = buffers[0].size();
    if (inverse)
    {
        leopard::ff8xor::IFFTButterflyBuffer(
            bytes, buffers[0].data(), buffers[1].data(), skew01);
        leopard::ff8xor::IFFTButterflyBuffer(
            bytes, buffers[2].data(), buffers[3].data(), skew23);
        leopard::ff8xor::IFFTButterflyBuffer(
            bytes, buffers[0].data(), buffers[2].data(), skew02);
        leopard::ff8xor::IFFTButterflyBuffer(
            bytes, buffers[1].data(), buffers[3].data(), skew02);
    }
    else
    {
        leopard::ff8xor::FFTButterflyBuffer(
            bytes, buffers[0].data(), buffers[2].data(), skew02);
        leopard::ff8xor::FFTButterflyBuffer(
            bytes, buffers[1].data(), buffers[3].data(), skew02);
        leopard::ff8xor::FFTButterflyBuffer(
            bytes, buffers[0].data(), buffers[1].data(), skew01);
        leopard::ff8xor::FFTButterflyBuffer(
            bytes, buffers[2].data(), buffers[3].data(), skew23);
    }
}

static bool TestAllGeneratedMaps()
{
    const uint64_t buffer_bytes = 1024;
    const leopard::ff8xor::FourBufferMode modes[] = {
        leopard::ff8xor::FourBufferMode::Xor2,
        leopard::ff8xor::FourBufferMode::Xor3
    };

    if (leopard::ff8xor::avx512four::GetTupleCount() != 64)
    {
        fprintf(stderr, "unexpected generated tuple count\n");
        return false;
    }

    for (unsigned tuple_index = 0; tuple_index < 64; ++tuple_index)
    {
        uint8_t skew01 = 0;
        uint8_t skew23 = 0;
        uint8_t skew02 = 0;
        if (!leopard::ff8xor::avx512four::GetTuple(
                tuple_index, skew01, skew23, skew02) ||
            leopard::ff8xor::avx512four::FindTupleIndex(
                skew01, skew23, skew02) != static_cast<int>(tuple_index))
        {
            fprintf(stderr, "tuple lookup mismatch: %u\n", tuple_index);
            return false;
        }

        FourBuffers input;
        for (unsigned index = 0; index < 4; ++index)
        {
            input[index].resize(static_cast<size_t>(buffer_bytes));
            Fill(input[index]);
        }

        for (unsigned direction = 0; direction < 2; ++direction)
        {
            const bool inverse = direction != 0;
            FourBuffers expected = input;
            leopard::ff8xor::SetFourBufferMode(
                leopard::ff8xor::FourBufferMode::Disabled);
            ApplyReference(expected, inverse, skew01, skew23, skew02);

            for (unsigned mode_index = 0; mode_index < 2; ++mode_index)
            {
                FourBuffers actual = input;
                leopard::ff8xor::SetFourBufferMode(modes[mode_index]);
                leopard::ff8xor::ResetFourBufferStatistics();
                if (!leopard::ff8xor::FourBufferButterflyBufferForTesting(
                        buffer_bytes,
                        actual[0].data(), actual[1].data(),
                        actual[2].data(), actual[3].data(),
                        inverse, skew01, skew23, skew02) ||
                    !EqualFour(
                        expected, actual, tuple_index, inverse,
                        modes[mode_index]))
                {
                    return false;
                }

                const leopard::ff8xor::FourBufferStatistics statistics =
                    leopard::ff8xor::GetLastFourBufferStatistics();
                const uint64_t original_factor =
                    (skew01 == 255 ? 3 : 4) +
                    (skew23 == 255 ? 3 : 4) +
                    (skew02 == 255 ? 6 : 8);
                if (statistics.FusedUnits != 1 ||
                    statistics.EstimatedPayloadBytesElided !=
                        buffer_bytes * (original_factor - 8))
                {
                    fprintf(stderr,
                        "four-buffer traffic statistics mismatch: tuple=%u\n",
                        tuple_index);
                    return false;
                }
            }
        }
    }

    return true;
}

static bool TestFallbacksDoNotTouchPayload()
{
    FourBuffers buffers;
    for (unsigned index = 0; index < 4; ++index)
    {
        buffers[index].resize(64);
        Fill(buffers[index]);
    }
    const FourBuffers original = buffers;
    leopard::ff8xor::SetFourBufferMode(
        leopard::ff8xor::FourBufferMode::Xor2);
    leopard::ff8xor::ResetFourBufferStatistics();
    if (leopard::ff8xor::FourBufferButterflyBufferForTesting(
            64,
            buffers[0].data(), buffers[1].data(),
            buffers[2].data(), buffers[3].data(),
            false, 1, 41, 26) || buffers != original)
    {
        fprintf(stderr, "narrow-tail fallback touched payload\n");
        return false;
    }
    FourBuffers partial;
    for (unsigned index = 0; index < 4; ++index)
    {
        partial[index].resize(960);
        Fill(partial[index]);
    }
    const FourBuffers partial_original = partial;
    if (leopard::ff8xor::FourBufferButterflyBufferForTesting(
            960,
            partial[0].data(), partial[1].data(),
            partial[2].data(), partial[3].data(),
            false, 1, 41, 26) || partial != partial_original)
    {
        fprintf(stderr, "partial-ZMM tail fallback touched payload\n");
        return false;
    }
    FourBuffers wide;
    for (unsigned index = 0; index < 4; ++index)
    {
        wide[index].resize(1024);
        Fill(wide[index]);
    }
    const FourBuffers wide_original = wide;
    if (leopard::ff8xor::FourBufferButterflyBufferForTesting(
            1024,
            wide[0].data(), wide[1].data(),
            wide[2].data(), wide[3].data(),
            false, 0, 0, 0) || wide != wide_original)
    {
        fprintf(stderr, "unknown-tuple fallback touched payload\n");
        return false;
    }
    leopard::ff8xor::SetFourBufferMode(
        leopard::ff8xor::FourBufferMode::Disabled);
    if (leopard::ff8xor::FourBufferButterflyBufferForTesting(
            1024,
            wide[0].data(), wide[1].data(),
            wide[2].data(), wide[3].data(),
            false, 1, 41, 26) || wide != wide_original)
    {
        fprintf(stderr, "disabled-mode fallback touched payload\n");
        return false;
    }
    const leopard::ff8xor::FourBufferStatistics statistics =
        leopard::ff8xor::GetLastFourBufferStatistics();
    if (statistics.FusedUnits != 0 ||
        statistics.EstimatedPayloadBytesElided != 0)
    {
        fprintf(stderr, "fallback changed fusion statistics\n");
        return false;
    }
    return true;
}

static bool TestNonZmmFallbackDoesNotTouchPayload()
{
    FourBuffers buffers;
    for (unsigned index = 0; index < 4; ++index)
    {
        buffers[index].resize(1024);
        Fill(buffers[index]);
    }
    const FourBuffers original = buffers;
    const leopard::ff8xor::KernelMode saved_kernel =
        leopard::ff8xor::GetKernelMode();
    const leopard::ff8xor::FourBufferMode saved_four =
        leopard::ff8xor::GetFourBufferMode();
    leopard::ff8xor::SetKernelMode(
        leopard::ff8xor::KernelMode::Portable);
    leopard::ff8xor::SetFourBufferMode(
        leopard::ff8xor::FourBufferMode::Xor3);
    leopard::ff8xor::ResetFourBufferStatistics();
    const bool applied =
        leopard::ff8xor::FourBufferButterflyBufferForTesting(
            1024,
            buffers[0].data(), buffers[1].data(),
            buffers[2].data(), buffers[3].data(),
            false, 1, 41, 26);
    const leopard::ff8xor::FourBufferStatistics statistics =
        leopard::ff8xor::GetLastFourBufferStatistics();
    leopard::ff8xor::SetFourBufferMode(saved_four);
    leopard::ff8xor::SetKernelMode(saved_kernel);

    if (applied || buffers != original || statistics.FusedUnits != 0 ||
        statistics.EstimatedPayloadBytesElided != 0)
    {
        fprintf(stderr, "non-ZMM fallback touched payload or statistics\n");
        return false;
    }
    return true;
}

struct Encoded
{
    unsigned OriginalCount;
    unsigned RecoveryCount;
    uint64_t BufferBytes;
    std::vector<Buffer> Original;
    std::vector<Buffer> Recovery;
};

static bool Encode(
    unsigned original_count,
    unsigned recovery_count,
    uint64_t buffer_bytes,
    leopard::ff8xor::FourBufferMode mode,
    const std::vector<Buffer>& original,
    std::vector<Buffer>& recovery)
{
    const unsigned work_count = leo_ff8xor_encode_work_count(
        original_count, recovery_count);
    std::vector<Buffer> work(work_count, Buffer(
        static_cast<size_t>(buffer_bytes), 0xa5));
    std::vector<const void*> original_ptrs(original_count);
    std::vector<void*> work_ptrs(work_count);
    for (unsigned i = 0; i < original_count; ++i)
        original_ptrs[i] = original[i].data();
    for (unsigned i = 0; i < work_count; ++i)
        work_ptrs[i] = work[i].data();

    leopard::ff8xor::SetFourBufferMode(mode);
    if (leo_ff8xor_encode(
            buffer_bytes, original_count, recovery_count, work_count,
            original_ptrs.data(), work_ptrs.data()) != Leopard_Success)
    {
        fprintf(stderr, "codec encode failed\n");
        return false;
    }
    recovery.assign(work.begin(), work.begin() + recovery_count);
    return true;
}

static bool Contains(
    const std::vector<unsigned>& indices,
    unsigned value)
{
    for (size_t i = 0; i < indices.size(); ++i)
    {
        if (indices[i] == value)
            return true;
    }
    return false;
}

static bool DecodePattern(
    const Encoded& encoded,
    leopard::ff8xor::FourBufferMode mode,
    const std::vector<unsigned>& missing_original,
    const std::vector<unsigned>& missing_recovery,
    bool require_fusion)
{
    const unsigned work_count = leo_ff8xor_decode_work_count(
        encoded.OriginalCount, encoded.RecoveryCount);
    std::vector<Buffer> work(work_count, Buffer(
        static_cast<size_t>(encoded.BufferBytes), 0x5a));
    std::vector<const void*> original_ptrs(encoded.OriginalCount);
    std::vector<const void*> recovery_ptrs(encoded.RecoveryCount);
    std::vector<void*> work_ptrs(work_count);
    for (unsigned i = 0; i < encoded.OriginalCount; ++i)
    {
        original_ptrs[i] = Contains(missing_original, i)
            ? NULL : encoded.Original[i].data();
    }
    for (unsigned i = 0; i < encoded.RecoveryCount; ++i)
    {
        recovery_ptrs[i] = Contains(missing_recovery, i)
            ? NULL : encoded.Recovery[i].data();
    }
    for (unsigned i = 0; i < work_count; ++i)
        work_ptrs[i] = work[i].data();

    leopard::ff8xor::SetFourBufferMode(mode);
    if (leo_ff8xor_decode(
            encoded.BufferBytes,
            encoded.OriginalCount,
            encoded.RecoveryCount,
            work_count,
            original_ptrs.data(),
            recovery_ptrs.data(),
            work_ptrs.data()) != Leopard_Success)
    {
        fprintf(stderr,
            "codec decode failed: k=%u r=%u bytes=%llu mode=%u\n",
            encoded.OriginalCount,
            encoded.RecoveryCount,
            static_cast<unsigned long long>(encoded.BufferBytes),
            static_cast<unsigned>(mode));
        return false;
    }
    for (size_t i = 0; i < missing_original.size(); ++i)
    {
        const unsigned missing = missing_original[i];
        if (work[missing] != encoded.Original[missing])
        {
            fprintf(stderr,
                "codec decode mismatch: missing=%u mode=%u\n",
                missing, static_cast<unsigned>(mode));
            return false;
        }
    }
    if (require_fusion &&
        leopard::ff8xor::GetLastFourBufferStatistics().FusedUnits == 0)
    {
        fprintf(stderr,
            "codec decode did not exercise fusion: k=%u r=%u bytes=%llu "
            "mode=%u losses=%u/%u\n",
            encoded.OriginalCount,
            encoded.RecoveryCount,
            static_cast<unsigned long long>(encoded.BufferBytes),
            static_cast<unsigned>(mode),
            static_cast<unsigned>(missing_original.size()),
            static_cast<unsigned>(missing_recovery.size()));
        return false;
    }
    return true;
}

static bool TestCodecEquivalence()
{
    const unsigned parameters[][3] = {
        { 4, 3, 512 },
        { 5, 3, 1024 },
        { 8, 4, 4096 },
        { 7, 3, 65536 },
        { 16, 4, 1024 },
        { 32, 8, 1024 },
        { 64, 16, 1024 },
        { 128, 32, 1024 },
        { 128, 128, 512 }
    };
    const leopard::ff8xor::FourBufferMode modes[] = {
        leopard::ff8xor::FourBufferMode::Xor2,
        leopard::ff8xor::FourBufferMode::Xor3
    };

    for (unsigned parameter_index = 0;
        parameter_index < sizeof(parameters) / sizeof(parameters[0]);
        ++parameter_index)
    {
        Encoded encoded;
        encoded.OriginalCount = parameters[parameter_index][0];
        encoded.RecoveryCount = parameters[parameter_index][1];
        encoded.BufferBytes = parameters[parameter_index][2];
        encoded.Original.resize(encoded.OriginalCount);
        for (unsigned i = 0; i < encoded.OriginalCount; ++i)
        {
            encoded.Original[i].resize(
                static_cast<size_t>(encoded.BufferBytes));
            Fill(encoded.Original[i]);
        }
        if (!Encode(
                encoded.OriginalCount,
                encoded.RecoveryCount,
                encoded.BufferBytes,
                leopard::ff8xor::FourBufferMode::Disabled,
                encoded.Original,
                encoded.Recovery))
        {
            return false;
        }

        for (unsigned mode_index = 0; mode_index < 2; ++mode_index)
        {
            std::vector<Buffer> recovery;
            if (!Encode(
                    encoded.OriginalCount,
                    encoded.RecoveryCount,
                    encoded.BufferBytes,
                    modes[mode_index],
                    encoded.Original,
                    recovery) || recovery != encoded.Recovery ||
                leopard::ff8xor::GetLastFourBufferStatistics().FusedUnits == 0)
            {
                fprintf(stderr,
                    "codec equivalence failed: parameter=%u mode=%u\n",
                    parameter_index,
                    static_cast<unsigned>(modes[mode_index]));
                return false;
            }

            std::vector<unsigned> one_original(
                1, encoded.OriginalCount / 2);
            const std::vector<unsigned> none;
            const bool require_decode_fusion =
                encoded.BufferBytes < 65536;
            if (!DecodePattern(
                    encoded, modes[mode_index], one_original, none,
                    require_decode_fusion))
                return false;

            std::vector<unsigned> maximum_original;
            for (unsigned i = 0; i < encoded.RecoveryCount; ++i)
                maximum_original.push_back(i);
            if (!DecodePattern(
                    encoded, modes[mode_index], maximum_original, none,
                    require_decode_fusion))
                return false;

            if (encoded.RecoveryCount >= 3)
            {
                std::vector<unsigned> mixed_original;
                mixed_original.push_back(0);
                mixed_original.push_back(encoded.OriginalCount - 1);
                std::vector<unsigned> mixed_recovery(1, 0);
                if (!DecodePattern(
                        encoded, modes[mode_index],
                        mixed_original, mixed_recovery,
                        require_decode_fusion))
                    return false;
            }
        }
    }
    return true;
}

static bool StatisticsAreZero()
{
    const leopard::ff8xor::FourBufferStatistics statistics =
        leopard::ff8xor::GetLastFourBufferStatistics();
    return statistics.FusedUnits == 0 &&
        statistics.EstimatedPayloadBytesElided == 0;
}

static bool TestPublicBoundaryStatisticsReset()
{
    const unsigned original_count = 4;
    const unsigned recovery_count = 3;
    const uint64_t buffer_bytes = 1024;
    std::vector<Buffer> original(original_count, Buffer(buffer_bytes));
    for (unsigned i = 0; i < original_count; ++i)
        Fill(original[i]);
    std::vector<Buffer> recovery;
    if (!Encode(
            original_count, recovery_count, buffer_bytes,
            leopard::ff8xor::FourBufferMode::Xor2,
            original, recovery) ||
        leopard::ff8xor::GetLastFourBufferStatistics().FusedUnits == 0)
    {
        fprintf(stderr, "failed to prime fusion statistics\n");
        return false;
    }

    if (leo_ff8xor_encode(1, 4, 3, 0, NULL, NULL) !=
            Leopard_InvalidSize || !StatisticsAreZero())
    {
        fprintf(stderr, "rejected encode retained fusion statistics\n");
        return false;
    }

    std::vector<const void*> one_original(1, original[0].data());
    std::vector<Buffer> one_work(1, Buffer(buffer_bytes));
    std::vector<void*> one_work_ptr(1, one_work[0].data());
    if (!Encode(
            original_count, recovery_count, buffer_bytes,
            leopard::ff8xor::FourBufferMode::Xor2,
            original, recovery) ||
        leo_ff8xor_encode(
            buffer_bytes, 1, 1, 1,
            one_original.data(), one_work_ptr.data()) != Leopard_Success ||
        !StatisticsAreZero())
    {
        fprintf(stderr, "k=1 encode retained fusion statistics\n");
        return false;
    }

    std::vector<const void*> originals(original_count);
    for (unsigned i = 0; i < original_count; ++i)
        originals[i] = original[i].data();
    const unsigned decode_count = leo_ff8xor_decode_work_count(
        original_count, recovery_count);
    std::vector<Buffer> decode_work(
        decode_count, Buffer(buffer_bytes));
    std::vector<void*> decode_ptrs(decode_count);
    std::vector<const void*> recovery_ptrs(recovery_count);
    for (unsigned i = 0; i < decode_count; ++i)
        decode_ptrs[i] = decode_work[i].data();
    if (!Encode(
            original_count, recovery_count, buffer_bytes,
            leopard::ff8xor::FourBufferMode::Xor2,
            original, recovery))
    {
        fprintf(stderr, "failed to re-prime no-loss decode statistics\n");
        return false;
    }
    // Encode replaces the Buffer objects, so form payload pointers only after
    // that replacement to avoid retaining invalidated vector storage.
    for (unsigned i = 0; i < recovery_count; ++i)
        recovery_ptrs[i] = recovery[i].data();
    if (leo_ff8xor_decode(
            buffer_bytes, original_count, recovery_count, decode_count,
            originals.data(), recovery_ptrs.data(), decode_ptrs.data()) !=
            Leopard_Success || !StatisticsAreZero())
    {
        fprintf(stderr, "no-loss decode retained fusion statistics\n");
        return false;
    }
    return true;
}

} // namespace

int main()
{
    if (leo_init() != 0)
    {
        fprintf(stderr, "leo_init failed\n");
        return 1;
    }
    // This fallback is required even in builds with no AVX-512 object.  Run
    // it before returning the hardware-specific CTest skip code.
    if (!TestNonZmmFallbackDoesNotTouchPayload())
        return 1;
    if (!leopard::ff8xor::IsKernelModeAvailable(
            leopard::ff8xor::KernelMode::Avx512Zmm))
    {
        printf("SKIP: AVX-512 ZMM kernels unavailable\n");
        return 77;
    }

    const leopard::ff8xor::KernelMode saved_kernel =
        leopard::ff8xor::GetKernelMode();
    const leopard::ff8xor::FourBufferMode saved_four =
        leopard::ff8xor::GetFourBufferMode();
    leopard::ff8xor::SetKernelMode(
        leopard::ff8xor::KernelMode::Avx512Zmm);

    const bool success = TestAllGeneratedMaps() &&
        TestFallbacksDoNotTouchPayload() &&
        TestCodecEquivalence() &&
        TestPublicBoundaryStatisticsReset();

    leopard::ff8xor::SetFourBufferMode(saved_four);
    leopard::ff8xor::SetKernelMode(saved_kernel);
    if (!success)
        return 1;
    printf("FF8 XOR AVX-512 four-buffer runtime tests passed\n");
    return 0;
}
