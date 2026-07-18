/*
    Finite deterministic tests for the fused FF8 XOR formal-derivative and
    top-sentinel-FFT boundary.
*/

#include "../LeopardFF8Xor.h"
#include "../leopard_ff8xor.h"

#include <algorithm>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <vector>

namespace {

typedef std::vector<uint8_t> Buffer;
typedef std::vector<Buffer> Buffers;

static uint32_t NextRandom(uint32_t& state)
{
    state ^= state << 13;
    state ^= state >> 17;
    state ^= state << 5;
    return state;
}

static void FillRandom(Buffers& buffers, uint32_t seed)
{
    for (size_t shard = 0; shard < buffers.size(); ++shard)
    {
        for (size_t offset = 0; offset < buffers[shard].size(); ++offset)
            buffers[shard][offset] = static_cast<uint8_t>(NextRandom(seed));
    }
}

static void Poison(Buffers& buffers, uint32_t seed)
{
    for (size_t shard = 0; shard < buffers.size(); ++shard)
    {
        for (size_t offset = 0; offset < buffers[shard].size(); ++offset)
        {
            uint8_t value = static_cast<uint8_t>(NextRandom(seed));
            if (value == 0)
                value = static_cast<uint8_t>(0x5aU ^ (shard & 0x3fU));
            buffers[shard][offset] = value;
        }
    }
}

static std::vector<void*> MutablePointers(Buffers& buffers)
{
    std::vector<void*> pointers(buffers.size());
    for (size_t index = 0; index < buffers.size(); ++index)
        pointers[index] = buffers[index].data();
    return pointers;
}

static std::vector<const void*> ConstPointers(const Buffers& buffers)
{
    std::vector<const void*> pointers(buffers.size());
    for (size_t index = 0; index < buffers.size(); ++index)
        pointers[index] = buffers[index].data();
    return pointers;
}

// This oracle deliberately preserves the old staged schedule instead of using
// the direct formula implemented by the candidate.  First apply the complete
// formal derivative, then the top skew-255 FFT layer (right ^= left).
static void ReferenceStagedBoundary(Buffers& work)
{
    const unsigned count = static_cast<unsigned>(work.size());
    const size_t buffer_bytes = work[0].size();
    for (unsigned index = 1; index < count; ++index)
    {
        const unsigned width = ((index ^ (index - 1)) + 1) >> 1;
        for (unsigned lane = 0; lane < width; ++lane)
        {
            Buffer& destination = work[index - width + lane];
            const Buffer& source = work[index + lane];
            for (size_t offset = 0; offset < buffer_bytes; ++offset)
                destination[offset] ^= source[offset];
        }
    }

    const unsigned half_count = count >> 1;
    for (unsigned index = 0; index < half_count; ++index)
    {
        for (size_t offset = 0; offset < buffer_bytes; ++offset)
            work[half_count + index][offset] ^= work[index][offset];
    }
}

static bool CheckBuffers(
    const Buffers& expected,
    const Buffers& actual,
    const char* label,
    leopard::ff8xor::KernelMode mode,
    unsigned count)
{
    if (expected == actual)
        return true;
    for (size_t shard = 0; shard < expected.size(); ++shard)
    {
        if (expected[shard] == actual[shard])
            continue;
        size_t offset = 0;
        while (offset < expected[shard].size() &&
               expected[shard][offset] == actual[shard][offset])
        {
            ++offset;
        }
        fprintf(stderr,
            "%s mismatch: mode=%u n=%u shard=%zu offset=%zu "
            "expected=%u actual=%u\n",
            label,
            static_cast<unsigned>(mode),
            count,
            shard,
            offset,
            static_cast<unsigned>(expected[shard][offset]),
            static_cast<unsigned>(actual[shard][offset]));
        return false;
    }
    return false;
}

static bool RunBoundaryCase(
    const Buffers& input,
    leopard::ff8xor::KernelMode mode,
    const char* label)
{
    Buffers expected = input;
    Buffers actual = input;
    ReferenceStagedBoundary(expected);
    std::vector<void*> pointers = MutablePointers(actual);
    leopard::ff8xor::FormalDerivativeTopFFTForTesting(
        actual[0].size(),
        static_cast<unsigned>(actual.size()),
        pointers.data());
    return CheckBuffers(
        expected, actual, label, mode, static_cast<unsigned>(actual.size()));
}

static bool TestExhaustiveLinearMaps(
    leopard::ff8xor::KernelMode mode)
{
    static const uint64_t kBufferBytes = 64;
    for (unsigned count = 2; count <= 256; count <<= 1)
    {
        Buffers pattern(count, Buffer(kBufferBytes, 0));
        if (!RunBoundaryCase(pattern, mode, "zero map"))
            return false;
        for (unsigned index = 0; index < count; ++index)
            std::fill(pattern[index].begin(), pattern[index].end(), 0xff);
        if (!RunBoundaryCase(pattern, mode, "all-one map"))
            return false;

        // XOR is bitwise and the map is binary-linear, so every one-hot shard
        // column proves the complete map for all bits and byte/vector lanes.
        for (unsigned basis = 0; basis < count; ++basis)
        {
            Buffers one_hot(count, Buffer(kBufferBytes, 0));
            std::fill(
                one_hot[basis].begin(), one_hot[basis].end(), 0xff);
            if (!RunBoundaryCase(one_hot, mode, "one-hot map"))
                return false;
        }

        uint32_t seed = UINT32_C(0xd3a10000) ^ count ^
            static_cast<unsigned>(mode);
        for (unsigned trial = 0; trial < 32; ++trial)
        {
            Buffers random(count, Buffer(kBufferBytes, 0));
            FillRandom(random, seed);
            if (!RunBoundaryCase(random, mode, "random map"))
                return false;
        }
    }
    return true;
}

static bool TestMixedWidthTailAndCanaries(
    leopard::ff8xor::KernelMode mode)
{
    static const unsigned kCount = 256;
    static const size_t kBufferBytes = 120;
    static const size_t kGuardBytes = 16;
    static const uint8_t kCanary = 0xcd;

    Buffers storage(
        kCount,
        Buffer(kBufferBytes + kGuardBytes * 2, kCanary));
    Buffers expected(kCount, Buffer(kBufferBytes, 0));
    uint32_t seed = UINT32_C(0x7a110000) ^ static_cast<unsigned>(mode);
    for (unsigned shard = 0; shard < kCount; ++shard)
    {
        for (size_t offset = 0; offset < kBufferBytes; ++offset)
        {
            const uint8_t value = static_cast<uint8_t>(NextRandom(seed));
            storage[shard][kGuardBytes + offset] = value;
            expected[shard][offset] = value;
        }
    }
    ReferenceStagedBoundary(expected);

    std::vector<void*> pointers(kCount);
    for (unsigned shard = 0; shard < kCount; ++shard)
        pointers[shard] = storage[shard].data() + kGuardBytes;
    leopard::ff8xor::FormalDerivativeTopFFTForTesting(
        kBufferBytes, kCount, pointers.data());

    Buffers actual(kCount, Buffer(kBufferBytes, 0));
    for (unsigned shard = 0; shard < kCount; ++shard)
    {
        memcpy(
            actual[shard].data(),
            storage[shard].data() + kGuardBytes,
            kBufferBytes);
        for (size_t guard = 0; guard < kGuardBytes; ++guard)
        {
            if (storage[shard][guard] != kCanary ||
                storage[shard][kGuardBytes + kBufferBytes + guard] != kCanary)
            {
                fprintf(stderr,
                    "boundary tail clobbered guard: mode=%u shard=%u\n",
                    static_cast<unsigned>(mode), shard);
                return false;
            }
        }
    }
    return CheckBuffers(expected, actual, "mixed-width tail", mode, kCount);
}

struct CodecCase
{
    unsigned OriginalCount;
    unsigned RecoveryCount;
    uint64_t BufferBytes;
    uint32_t Seed;
    const char* Label;
};

static bool RunDecode(
    const CodecCase& parameters,
    const Buffers& original,
    const Buffers& recovery,
    const std::vector<unsigned>& original_losses,
    const std::vector<unsigned>& recovery_losses)
{
    std::vector<const void*> original_ptrs = ConstPointers(original);
    std::vector<const void*> recovery_ptrs = ConstPointers(recovery);
    for (size_t index = 0; index < original_losses.size(); ++index)
        original_ptrs[original_losses[index]] = NULL;
    for (size_t index = 0; index < recovery_losses.size(); ++index)
        recovery_ptrs[recovery_losses[index]] = NULL;

    const unsigned work_count = leo_ff8xor_decode_work_count(
        parameters.OriginalCount, parameters.RecoveryCount);
    Buffers work(
        work_count,
        Buffer(static_cast<size_t>(parameters.BufferBytes), 0));
    Poison(work, parameters.Seed ^ UINT32_C(0xdec0de00));
    std::vector<void*> work_ptrs = MutablePointers(work);
    const LeopardResult result = leo_ff8xor_decode(
        parameters.BufferBytes,
        parameters.OriginalCount,
        parameters.RecoveryCount,
        work_count,
        original_ptrs.data(),
        recovery_ptrs.data(),
        work_ptrs.data());
    if (result != Leopard_Success)
    {
        fprintf(stderr, "%s decode failed: %s\n",
            parameters.Label, leo_result_string(result));
        return false;
    }
    for (size_t index = 0; index < original_losses.size(); ++index)
    {
        const unsigned shard = original_losses[index];
        if (work[shard] != original[shard])
        {
            fprintf(stderr, "%s recovery mismatch: shard=%u\n",
                parameters.Label, shard);
            return false;
        }
    }
    return true;
}

static bool RunCodecCase(const CodecCase& parameters)
{
    Buffers original(
        parameters.OriginalCount,
        Buffer(static_cast<size_t>(parameters.BufferBytes), 0));
    FillRandom(original, parameters.Seed);
    const unsigned encode_work_count = leo_ff8xor_encode_work_count(
        parameters.OriginalCount, parameters.RecoveryCount);
    Buffers encode_work(
        encode_work_count,
        Buffer(static_cast<size_t>(parameters.BufferBytes), 0));
    Poison(encode_work, parameters.Seed ^ UINT32_C(0xface0000));
    std::vector<const void*> original_ptrs = ConstPointers(original);
    std::vector<void*> encode_ptrs = MutablePointers(encode_work);
    const LeopardResult encode_result = leo_ff8xor_encode(
        parameters.BufferBytes,
        parameters.OriginalCount,
        parameters.RecoveryCount,
        encode_work_count,
        original_ptrs.data(),
        encode_ptrs.data());
    if (encode_result != Leopard_Success)
    {
        fprintf(stderr, "%s encode failed: %s\n",
            parameters.Label, leo_result_string(encode_result));
        return false;
    }
    Buffers recovery(
        encode_work.begin(),
        encode_work.begin() + parameters.RecoveryCount);

    std::vector<unsigned> sparse_original(1, parameters.OriginalCount - 1);
    std::vector<unsigned> sparse_recovery;
    for (unsigned index = 1; index < parameters.RecoveryCount; ++index)
        sparse_recovery.push_back(index);
    if (!RunDecode(
            parameters,
            original,
            recovery,
            sparse_original,
            sparse_recovery))
    {
        return false;
    }

    const unsigned mixed_count = parameters.RecoveryCount < 3
        ? parameters.RecoveryCount
        : 3;
    std::vector<unsigned> mixed_original(mixed_count);
    for (unsigned index = 0; index < mixed_count; ++index)
        mixed_original[index] = (index * 7 + 1) % parameters.OriginalCount;
    std::sort(mixed_original.begin(), mixed_original.end());
    mixed_original.erase(
        std::unique(mixed_original.begin(), mixed_original.end()),
        mixed_original.end());
    std::vector<unsigned> mixed_recovery;
    for (unsigned index = mixed_count; index < parameters.RecoveryCount; ++index)
        mixed_recovery.push_back(index);
    if (!RunDecode(
            parameters,
            original,
            recovery,
            mixed_original,
            mixed_recovery))
    {
        return false;
    }

    std::vector<unsigned> maximum_original(parameters.RecoveryCount);
    for (unsigned index = 0; index < parameters.RecoveryCount; ++index)
        maximum_original[index] = index;
    const std::vector<unsigned> no_recovery_losses;
    return RunDecode(
        parameters,
        original,
        recovery,
        maximum_original,
        no_recovery_losses);
}

static bool TestEndToEnd(leopard::ff8xor::KernelMode mode)
{
    // Full and truncated transforms cover both untracked and tracked decoder
    // paths.  Sparse, mixed, and maximum losses drive the ErrorBitfield at
    // every remaining mip after the fused top layer.
    static const CodecCase kCases[] = {
        { 2, 2, 64, UINT32_C(0x04020002), "full n=4" },
        { 9, 3, 64, UINT32_C(0x10090003), "truncated n=16" },
        { 8, 8, 64, UINT32_C(0x10080008), "full n=16" },
        { 128, 128, 64, UINT32_C(0x80128080), "full n=256" },
        { 9, 3, 65536, UINT32_C(0x64090003),
            "tracked poisoned truncated n=16" }
    };
    for (size_t index = 0; index < sizeof(kCases) / sizeof(kCases[0]); ++index)
    {
        CodecCase parameters = kCases[index];
        parameters.Seed ^= static_cast<unsigned>(mode);
        if (!RunCodecCase(parameters))
            return false;
    }
    return true;
}

static bool TestEveryPaddedSize()
{
    static const CodecCase kCases[] = {
        { 2, 2, 64, UINT32_C(0x04000001), "padded n=4" },
        { 5, 2, 64, UINT32_C(0x08000001), "padded n=8" },
        { 9, 3, 64, UINT32_C(0x10000001), "padded n=16" },
        { 17, 5, 64, UINT32_C(0x20000001), "padded n=32" },
        { 33, 9, 64, UINT32_C(0x40000001), "padded n=64" },
        { 65, 17, 64, UINT32_C(0x80000001), "padded n=128" },
        { 128, 64, 64, UINT32_C(0x25600001), "padded n=256" },
        // A small recovery power with k=128 gives a highly truncated n=256
        // transform and many poisoned deferred-zero buffers without making the
        // sanitizer case needlessly allocate 128 available recovery shards.
        { 128, 2, 65536, UINT32_C(0x25664002),
            "tracked poisoned truncated n=256" }
    };
    for (size_t index = 0; index < sizeof(kCases) / sizeof(kCases[0]); ++index)
    {
        if (!RunCodecCase(kCases[index]))
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

    typedef leopard::ff8xor::KernelMode KernelMode;
    static const KernelMode kModes[] = {
        KernelMode::Portable,
        KernelMode::Simd128,
        KernelMode::Avx2,
        KernelMode::Avx512VL,
        KernelMode::Avx512Zmm
    };
    const KernelMode saved_mode = leopard::ff8xor::GetKernelMode();
    const leopard::ff8xor::FourBufferMode saved_four =
        leopard::ff8xor::GetFourBufferMode();
    leopard::ff8xor::SetFourBufferMode(
        leopard::ff8xor::FourBufferMode::Disabled);

    bool success = true;
    for (size_t index = 0; index < sizeof(kModes) / sizeof(kModes[0]); ++index)
    {
        const KernelMode mode = kModes[index];
        if (!leopard::ff8xor::IsKernelModeAvailable(mode))
            continue;
        leopard::ff8xor::SetKernelMode(mode);
        if (leopard::ff8xor::GetActiveKernelMode() != mode ||
            !TestExhaustiveLinearMaps(mode) ||
            !TestMixedWidthTailAndCanaries(mode) ||
            !TestEndToEnd(mode))
        {
            success = false;
            break;
        }
    }

    if (success)
    {
        leopard::ff8xor::SetKernelMode(KernelMode::Auto);
        success = TestEveryPaddedSize();
    }

    leopard::ff8xor::SetFourBufferMode(saved_four);
    leopard::ff8xor::SetKernelMode(saved_mode);
    if (!success)
        return 1;
    printf("FF8 XOR formal-derivative boundary tests passed\n");
    return 0;
}
