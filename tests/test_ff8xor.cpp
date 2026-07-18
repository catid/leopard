/*
    Finite deterministic tests for the experimental FF8 XOR-circuit backend.
*/

#include "../LeopardCommon.h"
#include "../LeopardFF8Xor.h"
#include "../leopard_ff8xor.h"

#include <algorithm>
#include <limits.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <vector>

namespace leopard { namespace ff8xor {

static LEO_FORCE_INLINE uint8_t XorValue(uint8_t a, uint8_t b)
{
    return static_cast<uint8_t>(a ^ b);
}

}} // namespace leopard::ff8xor

#include "../generated/LeopardFF8XorCircuits.inl"

namespace {

static const unsigned kFieldBits = 8;
static const unsigned kFieldOrder = 256;
static const unsigned kFieldModulus = 255;
static const unsigned kFieldPolynomial = 0x11d;
static const uint8_t kCantorBasis[kFieldBits] = {
    1, 214, 152, 146, 86, 200, 88, 230
};

static uint8_t LogLUT[kFieldOrder];
static uint8_t ExpLUT[kFieldOrder];

static uint8_t AddMod(uint8_t a, uint8_t b)
{
    const unsigned sum = static_cast<unsigned>(a) + b;
    return static_cast<uint8_t>(sum + (sum >> kFieldBits));
}

static void InitializeReferenceField()
{
    unsigned state = 1;
    for (unsigned exponent = 0; exponent < kFieldModulus; ++exponent)
    {
        ExpLUT[state] = static_cast<uint8_t>(exponent);
        state <<= 1;
        if (state >= kFieldOrder)
            state ^= kFieldPolynomial;
    }
    ExpLUT[0] = static_cast<uint8_t>(kFieldModulus);

    LogLUT[0] = 0;
    for (unsigned bit = 0; bit < kFieldBits; ++bit)
    {
        const unsigned width = 1U << bit;
        for (unsigned value = 0; value < width; ++value)
            LogLUT[value + width] = LogLUT[value] ^ kCantorBasis[bit];
    }
    for (unsigned value = 0; value < kFieldOrder; ++value)
        LogLUT[value] = ExpLUT[LogLUT[value]];
    for (unsigned value = 0; value < kFieldOrder; ++value)
        ExpLUT[LogLUT[value]] = static_cast<uint8_t>(value);
    ExpLUT[kFieldModulus] = ExpLUT[0];
}

static uint8_t ReferenceMultiply(uint8_t value, uint8_t log_multiplier)
{
    if (value == 0)
        return 0;
    return ExpLUT[AddMod(LogLUT[value], log_multiplier)];
}

static uint16_t ReferenceFFT(uint16_t state, unsigned skew)
{
    uint8_t x = static_cast<uint8_t>(state);
    uint8_t y = static_cast<uint8_t>(state >> 8);
    if (skew != kFieldModulus)
        x ^= ReferenceMultiply(y, static_cast<uint8_t>(skew));
    y ^= x;
    return static_cast<uint16_t>(x | (static_cast<uint16_t>(y) << 8));
}

static uint16_t ReferenceIFFT(uint16_t state, unsigned skew)
{
    uint8_t x = static_cast<uint8_t>(state);
    uint8_t y = static_cast<uint8_t>(state >> 8);
    y ^= x;
    if (skew != kFieldModulus)
        x ^= ReferenceMultiply(y, static_cast<uint8_t>(skew));
    return static_cast<uint16_t>(x | (static_cast<uint16_t>(y) << 8));
}

template <unsigned Coefficient>
static uint8_t ApplyMultiplyCircuit(uint8_t value)
{
    uint8_t x0 = static_cast<uint8_t>((value >> 0) & 1);
    uint8_t x1 = static_cast<uint8_t>((value >> 1) & 1);
    uint8_t x2 = static_cast<uint8_t>((value >> 2) & 1);
    uint8_t x3 = static_cast<uint8_t>((value >> 3) & 1);
    uint8_t x4 = static_cast<uint8_t>((value >> 4) & 1);
    uint8_t x5 = static_cast<uint8_t>((value >> 5) & 1);
    uint8_t x6 = static_cast<uint8_t>((value >> 6) & 1);
    uint8_t x7 = static_cast<uint8_t>((value >> 7) & 1);

    leopard::ff8xor::generated::MultiplyCircuit<Coefficient>::Apply(
        x0, x1, x2, x3, x4, x5, x6, x7);

    return static_cast<uint8_t>(
        (x0 << 0) | (x1 << 1) | (x2 << 2) | (x3 << 3) |
        (x4 << 4) | (x5 << 5) | (x6 << 6) | (x7 << 7));
}

template <unsigned Skew, bool Inverse>
static uint16_t ApplyButterflyCircuit(uint16_t state)
{
    uint8_t x0 = static_cast<uint8_t>((state >> 0) & 1);
    uint8_t x1 = static_cast<uint8_t>((state >> 1) & 1);
    uint8_t x2 = static_cast<uint8_t>((state >> 2) & 1);
    uint8_t x3 = static_cast<uint8_t>((state >> 3) & 1);
    uint8_t x4 = static_cast<uint8_t>((state >> 4) & 1);
    uint8_t x5 = static_cast<uint8_t>((state >> 5) & 1);
    uint8_t x6 = static_cast<uint8_t>((state >> 6) & 1);
    uint8_t x7 = static_cast<uint8_t>((state >> 7) & 1);
    uint8_t y0 = static_cast<uint8_t>((state >> 8) & 1);
    uint8_t y1 = static_cast<uint8_t>((state >> 9) & 1);
    uint8_t y2 = static_cast<uint8_t>((state >> 10) & 1);
    uint8_t y3 = static_cast<uint8_t>((state >> 11) & 1);
    uint8_t y4 = static_cast<uint8_t>((state >> 12) & 1);
    uint8_t y5 = static_cast<uint8_t>((state >> 13) & 1);
    uint8_t y6 = static_cast<uint8_t>((state >> 14) & 1);
    uint8_t y7 = static_cast<uint8_t>((state >> 15) & 1);

    if (Inverse)
    {
        leopard::ff8xor::generated::IFFTCircuit<Skew>::Apply(
            x0, x1, x2, x3, x4, x5, x6, x7,
            y0, y1, y2, y3, y4, y5, y6, y7);
    }
    else
    {
        leopard::ff8xor::generated::FFTCircuit<Skew>::Apply(
            x0, x1, x2, x3, x4, x5, x6, x7,
            y0, y1, y2, y3, y4, y5, y6, y7);
    }

    return static_cast<uint16_t>(
        (x0 << 0) | (x1 << 1) | (x2 << 2) | (x3 << 3) |
        (x4 << 4) | (x5 << 5) | (x6 << 6) | (x7 << 7) |
        (static_cast<uint16_t>(y0) << 8) |
        (static_cast<uint16_t>(y1) << 9) |
        (static_cast<uint16_t>(y2) << 10) |
        (static_cast<uint16_t>(y3) << 11) |
        (static_cast<uint16_t>(y4) << 12) |
        (static_cast<uint16_t>(y5) << 13) |
        (static_cast<uint16_t>(y6) << 14) |
        (static_cast<uint16_t>(y7) << 15));
}

typedef uint8_t (*MultiplyCircuitFunction)(uint8_t);
typedef uint16_t (*ButterflyCircuitFunction)(uint16_t);

#define LEO_MULTIPLY_ENTRY(Coefficient) &ApplyMultiplyCircuit<Coefficient>,
static const MultiplyCircuitFunction MultiplyCircuits[kFieldOrder] = {
    LEOPARD_FF8XOR_FOR_EACH_LOG(LEO_MULTIPLY_ENTRY)
};
#undef LEO_MULTIPLY_ENTRY

#define LEO_FFT_ENTRY(Skew) &ApplyButterflyCircuit<Skew, false>,
static const ButterflyCircuitFunction FFTCircuits[kFieldOrder] = {
    LEOPARD_FF8XOR_FOR_EACH_LOG(LEO_FFT_ENTRY)
};
#undef LEO_FFT_ENTRY

#define LEO_IFFT_ENTRY(Skew) &ApplyButterflyCircuit<Skew, true>,
static const ButterflyCircuitFunction IFFTCircuits[kFieldOrder] = {
    LEOPARD_FF8XOR_FOR_EACH_LOG(LEO_IFFT_ENTRY)
};
#undef LEO_IFFT_ENTRY

static bool TestMultiplicationCircuits()
{
    for (unsigned log_multiplier = 0; log_multiplier < kFieldOrder; ++log_multiplier)
    {
        for (unsigned value = 0; value < kFieldOrder; ++value)
        {
            const uint8_t expected = ReferenceMultiply(
                static_cast<uint8_t>(value),
                static_cast<uint8_t>(log_multiplier));
            const uint8_t actual = MultiplyCircuits[log_multiplier](
                static_cast<uint8_t>(value));
            if (actual != expected)
            {
                fprintf(stderr,
                    "multiply mismatch: log=%u value=%u expected=%u actual=%u\n",
                    log_multiplier, value, expected, actual);
                return false;
            }
        }
    }
    return true;
}

static uint32_t NextRandom(uint32_t& state)
{
    state ^= state << 13;
    state ^= state >> 17;
    state ^= state << 5;
    return state;
}

static bool CheckButterflyState(unsigned skew, uint16_t state)
{
    const uint16_t expected_fft = ReferenceFFT(state, skew);
    const uint16_t expected_ifft = ReferenceIFFT(state, skew);
    const uint16_t actual_fft = FFTCircuits[skew](state);
    const uint16_t actual_ifft = IFFTCircuits[skew](state);

    if (actual_fft != expected_fft || actual_ifft != expected_ifft)
    {
        fprintf(stderr,
            "butterfly mismatch: skew=%u state=%04x fft=%04x/%04x ifft=%04x/%04x\n",
            skew, state, actual_fft, expected_fft, actual_ifft, expected_ifft);
        return false;
    }
    if (IFFTCircuits[skew](actual_fft) != state ||
        FFTCircuits[skew](actual_ifft) != state)
    {
        fprintf(stderr,
            "butterfly inverse mismatch: skew=%u state=%04x\n",
            skew, state);
        return false;
    }
    return true;
}

static bool TestButterflyCircuits()
{
    uint32_t random_state = 0xff8c1ac0U;
    for (unsigned skew = 0; skew < kFieldOrder; ++skew)
    {
        if (!CheckButterflyState(skew, 0) ||
            !CheckButterflyState(skew, 0xffff))
            return false;

        for (unsigned bit = 0; bit < 16; ++bit)
            if (!CheckButterflyState(skew, static_cast<uint16_t>(1U << bit)))
                return false;

        for (unsigned trial = 0; trial < 128; ++trial)
        {
            const uint16_t state = static_cast<uint16_t>(NextRandom(random_state));
            if (!CheckButterflyState(skew, state))
                return false;
        }
    }
    return true;
}


//------------------------------------------------------------------------------
// Packed/plane helpers

typedef std::vector<uint8_t> Buffer;
typedef std::vector<Buffer> Buffers;

static void PackedToPlane(const Buffer& packed, Buffer& plane)
{
    const size_t buffer_bytes = packed.size();
    const size_t plane_bytes = buffer_bytes / kFieldBits;
    plane.assign(buffer_bytes, 0);

    for (size_t group = 0; group < plane_bytes; ++group)
    {
        for (unsigned lane = 0; lane < kFieldBits; ++lane)
        {
            const uint8_t packed_value = packed[group * kFieldBits + lane];
            for (unsigned bit = 0; bit < kFieldBits; ++bit)
            {
                plane[bit * plane_bytes + group] |= static_cast<uint8_t>(
                    ((packed_value >> bit) & 1U) << lane);
            }
        }
    }
}

static void PlaneToPacked(const Buffer& plane, Buffer& packed)
{
    const size_t buffer_bytes = plane.size();
    const size_t plane_bytes = buffer_bytes / kFieldBits;
    packed.assign(buffer_bytes, 0);

    for (size_t group = 0; group < plane_bytes; ++group)
    {
        for (unsigned lane = 0; lane < kFieldBits; ++lane)
        {
            uint8_t packed_value = 0;
            for (unsigned bit = 0; bit < kFieldBits; ++bit)
            {
                packed_value |= static_cast<uint8_t>(
                    ((plane[bit * plane_bytes + group] >> lane) & 1U) << bit);
            }
            packed[group * kFieldBits + lane] = packed_value;
        }
    }
}

static void FillRandom(Buffer& buffer, uint32_t& state)
{
    for (size_t index = 0; index < buffer.size(); ++index)
        buffer[index] = static_cast<uint8_t>(NextRandom(state));
}

static void FillRandom(Buffers& buffers, uint32_t seed)
{
    for (size_t index = 0; index < buffers.size(); ++index)
        FillRandom(buffers[index], seed);
}

static Buffers AllocateBuffers(unsigned count, uint64_t buffer_bytes)
{
    return Buffers(count, Buffer(static_cast<size_t>(buffer_bytes), 0));
}

static std::vector<const void*> GetConstPointers(const Buffers& buffers)
{
    std::vector<const void*> pointers(buffers.size());
    for (size_t index = 0; index < buffers.size(); ++index)
        pointers[index] = buffers[index].data();
    return pointers;
}

static std::vector<void*> GetMutablePointers(Buffers& buffers)
{
    std::vector<void*> pointers(buffers.size());
    for (size_t index = 0; index < buffers.size(); ++index)
        pointers[index] = buffers[index].data();
    return pointers;
}

static bool CheckEqual(
    const Buffer& actual,
    const Buffer& expected,
    const char* operation,
    unsigned index)
{
    if (actual.size() == expected.size() &&
        memcmp(actual.data(), expected.data(), actual.size()) == 0)
        return true;

    size_t mismatch = 0;
    const size_t common = actual.size() < expected.size()
        ? actual.size()
        : expected.size();
    while (mismatch < common && actual[mismatch] == expected[mismatch])
        ++mismatch;

    fprintf(stderr,
        "%s buffer %u mismatch at byte %llu (sizes %llu/%llu)\n",
        operation,
        index,
        static_cast<unsigned long long>(mismatch),
        static_cast<unsigned long long>(actual.size()),
        static_cast<unsigned long long>(expected.size()));
    return false;
}

static bool TestTransposeHelpers()
{
    static const uint64_t kSizes[] = { 64, 192, 1024, 4096 };
    uint32_t state = 0x8badf00dU;

    for (unsigned size_index = 0;
         size_index < sizeof(kSizes) / sizeof(kSizes[0]);
         ++size_index)
    {
        Buffer packed(static_cast<size_t>(kSizes[size_index]));
        Buffer plane;
        Buffer round_trip;
        FillRandom(packed, state);
        PackedToPlane(packed, plane);

        const size_t plane_bytes = plane.size() / kFieldBits;
        for (size_t group = 0; group < plane_bytes; ++group)
        {
            for (unsigned lane = 0; lane < kFieldBits; ++lane)
            {
                for (unsigned bit = 0; bit < kFieldBits; ++bit)
                {
                    const unsigned packed_bit =
                        (packed[group * kFieldBits + lane] >> bit) & 1U;
                    const unsigned plane_bit =
                        (plane[bit * plane_bytes + group] >> lane) & 1U;
                    if (packed_bit != plane_bit)
                    {
                        fprintf(stderr,
                            "transpose formula mismatch: bytes=%llu group=%llu lane=%u bit=%u\n",
                            static_cast<unsigned long long>(kSizes[size_index]),
                            static_cast<unsigned long long>(group),
                            lane,
                            bit);
                        return false;
                    }
                }
            }
        }

        PlaneToPacked(plane, round_trip);
        if (!CheckEqual(round_trip, packed, "transpose round trip", size_index))
            return false;
    }
    return true;
}


//------------------------------------------------------------------------------
// Public API validation

static bool ExpectResult(
    LeopardResult actual,
    LeopardResult expected,
    const char* operation)
{
    if (actual == expected)
        return true;
    fprintf(stderr,
        "%s returned %d (%s), expected %d (%s)\n",
        operation,
        static_cast<int>(actual),
        leo_result_string(actual),
        static_cast<int>(expected),
        leo_result_string(expected));
    return false;
}

static bool TestWorkCounts()
{
    static const unsigned kCounts[][2] = {
        { 1, 1 }, { 8, 1 }, { 4, 2 }, { 16, 4 },
        { 64, 16 }, { 128, 32 }, { 128, 128 }
    };
    for (unsigned index = 0;
         index < sizeof(kCounts) / sizeof(kCounts[0]);
         ++index)
    {
        const unsigned original_count = kCounts[index][0];
        const unsigned recovery_count = kCounts[index][1];
        if (leo_ff8xor_encode_work_count(original_count, recovery_count) !=
                leo_encode_work_count(original_count, recovery_count) ||
            leo_ff8xor_decode_work_count(original_count, recovery_count) !=
                leo_decode_work_count(original_count, recovery_count))
        {
            fprintf(stderr,
                "work-count mismatch: k=%u r=%u\n",
                original_count,
                recovery_count);
            return false;
        }
    }

    static const unsigned kInvalidCounts[][2] = {
        { 0, 0 }, { 4, 0 }, { 4, 5 }, { 129, 128 }, { 256, 1 }
    };
    for (unsigned index = 0;
         index < sizeof(kInvalidCounts) / sizeof(kInvalidCounts[0]);
         ++index)
    {
        const unsigned original_count = kInvalidCounts[index][0];
        const unsigned recovery_count = kInvalidCounts[index][1];
        if (leo_ff8xor_encode_work_count(
                original_count, recovery_count) != 0 ||
            leo_ff8xor_decode_work_count(
                original_count, recovery_count) != 0)
        {
            fprintf(stderr,
                "invalid work count was nonzero: k=%u r=%u\n",
                original_count,
                recovery_count);
            return false;
        }
    }
    return true;
}

static bool TestCallInitialize()
{
    Buffers original = AllocateBuffers(2, 64);
    Buffers recovery = AllocateBuffers(2, 64);
    Buffers work = AllocateBuffers(4, 64);
    std::vector<const void*> original_ptrs = GetConstPointers(original);
    std::vector<const void*> recovery_ptrs = GetConstPointers(recovery);
    std::vector<void*> work_ptrs = GetMutablePointers(work);

    return ExpectResult(
               leo_ff8xor_encode(
                   64, 2, 2, 4, original_ptrs.data(), work_ptrs.data()),
               Leopard_CallInitialize,
               "encode before leo_init") &&
           ExpectResult(
               leo_ff8xor_decode(
                   64, 2, 2, 4,
                   original_ptrs.data(), recovery_ptrs.data(), work_ptrs.data()),
               Leopard_CallInitialize,
               "decode before leo_init");
}

static bool TestAPIValidation()
{
    Buffers original = AllocateBuffers(4, 64);
    Buffers recovery = AllocateBuffers(2, 64);
    Buffers work = AllocateBuffers(8, 64);
    std::vector<const void*> original_ptrs = GetConstPointers(original);
    std::vector<const void*> recovery_ptrs = GetConstPointers(recovery);
    std::vector<void*> work_ptrs = GetMutablePointers(work);

    if (!ExpectResult(
            leo_ff8xor_encode(
                0, 4, 2, 4, original_ptrs.data(), work_ptrs.data()),
            Leopard_InvalidSize,
            "encode zero bytes") ||
        !ExpectResult(
            leo_ff8xor_encode(
                63, 4, 2, 4, original_ptrs.data(), work_ptrs.data()),
            Leopard_InvalidSize,
            "encode non-multiple size") ||
        !ExpectResult(
            leo_ff8xor_encode(
                64, 4, 0, 4, original_ptrs.data(), work_ptrs.data()),
            Leopard_InvalidCounts,
            "encode zero recovery") ||
        !ExpectResult(
            leo_ff8xor_encode(
                64, 4, 5, 4, original_ptrs.data(), work_ptrs.data()),
            Leopard_InvalidCounts,
            "encode recovery greater than original") ||
        !ExpectResult(
            leo_ff8xor_encode(64, 4, 2, 4, NULL, work_ptrs.data()),
            Leopard_InvalidInput,
            "encode null originals") ||
        !ExpectResult(
            leo_ff8xor_encode(64, 4, 2, 4, original_ptrs.data(), NULL),
            Leopard_InvalidInput,
            "encode null work") ||
        !ExpectResult(
            leo_ff8xor_encode(
                64, 4, 2, 3, original_ptrs.data(), work_ptrs.data()),
            Leopard_InvalidCounts,
            "encode wrong work count") ||
        !ExpectResult(
            leo_ff8xor_encode(
                64, 4, 1, 0, original_ptrs.data(), work_ptrs.data()),
            Leopard_InvalidCounts,
            "encode r=1 wrong work count") ||
        !ExpectResult(
            leo_ff8xor_encode(
                64, 1, 1, 0, original_ptrs.data(), work_ptrs.data()),
            Leopard_InvalidCounts,
            "encode k=1 wrong work count"))
        return false;

    if (!ExpectResult(
            leo_ff8xor_decode(
                0, 4, 2, 8,
                original_ptrs.data(), recovery_ptrs.data(), work_ptrs.data()),
            Leopard_InvalidSize,
            "decode zero bytes") ||
        !ExpectResult(
            leo_ff8xor_decode(
                64, 4, 0, 8,
                original_ptrs.data(), recovery_ptrs.data(), work_ptrs.data()),
            Leopard_InvalidCounts,
            "decode zero recovery") ||
        !ExpectResult(
            leo_ff8xor_decode(
                64, 4, 2, 8,
                NULL, recovery_ptrs.data(), work_ptrs.data()),
            Leopard_InvalidInput,
            "decode null originals") ||
        !ExpectResult(
            leo_ff8xor_decode(
                64, 4, 2, 8,
                original_ptrs.data(), NULL, work_ptrs.data()),
            Leopard_InvalidInput,
            "decode null recovery") ||
        !ExpectResult(
            leo_ff8xor_decode(
                64, 4, 2, 8,
                original_ptrs.data(), recovery_ptrs.data(), NULL),
            Leopard_InvalidInput,
            "decode null work") ||
        !ExpectResult(
            leo_ff8xor_decode(
                64, 4, 2, 7,
                original_ptrs.data(), recovery_ptrs.data(), work_ptrs.data()),
            Leopard_InvalidCounts,
            "decode no-loss wrong work count"))
        return false;

    // Force the general decoder path so the work count is checked.
    original_ptrs[0] = NULL;
    if (!ExpectResult(
            leo_ff8xor_decode(
                64, 4, 2, 7,
                original_ptrs.data(), recovery_ptrs.data(), work_ptrs.data()),
            Leopard_InvalidCounts,
            "decode wrong work count"))
        return false;
    original_ptrs[0] = original[0].data();

    // A no-loss decode does not require any recovery shard.  In particular,
    // k=1 must copy the available original instead of dereferencing a missing
    // recovery pointer.
    Buffers singleton_original = AllocateBuffers(1, 64);
    Buffers singleton_work = AllocateBuffers(1, 64);
    uint32_t singleton_seed = 0x1dec0deU;
    FillRandom(singleton_original[0], singleton_seed);
    std::vector<const void*> singleton_original_ptrs =
        GetConstPointers(singleton_original);
    std::vector<const void*> singleton_recovery_ptrs(1, NULL);
    std::vector<void*> singleton_work_ptrs =
        GetMutablePointers(singleton_work);
    if (!ExpectResult(
            leo_ff8xor_decode(
                64, 1, 1, 0,
                singleton_original_ptrs.data(),
                singleton_recovery_ptrs.data(),
                singleton_work_ptrs.data()),
            Leopard_InvalidCounts,
            "decode k=1 wrong work count") ||
        !ExpectResult(
            leo_ff8xor_decode(
                64, 1, 1, 1,
                singleton_original_ptrs.data(),
                singleton_recovery_ptrs.data(),
                singleton_work_ptrs.data()),
            Leopard_Success,
            "decode k=1 no loss with missing recovery") ||
        !CheckEqual(
            singleton_work[0], singleton_original[0],
            "decode k=1 no loss with missing recovery", 0))
        return false;

    Buffer dummy(64, 0);
    std::vector<const void*> large_original(129, dummy.data());
    std::vector<const void*> large_recovery(128, dummy.data());
    // Unsupported helper queries return zero.  Supply one non-null dummy slot
    // here so the operation reaches and reports the FF8 range check.
    std::vector<void*> large_encode_work(1, dummy.data());
    std::vector<void*> large_decode_work(1, dummy.data());
    const void* bounded_original[] = { dummy.data() };
    const void* bounded_recovery[] = { dummy.data() };
    void* bounded_work[] = { dummy.data() };

    return ExpectResult(
               leo_ff8xor_encode(
                   64, 129, 128,
                   static_cast<unsigned>(large_encode_work.size()),
                   large_original.data(), large_encode_work.data()),
               Leopard_TooMuchData,
               "encode n greater than 256") &&
           ExpectResult(
               leo_ff8xor_decode(
                   64, 129, 128,
                   static_cast<unsigned>(large_decode_work.size()),
                   large_original.data(), large_recovery.data(),
                   large_decode_work.data()),
               Leopard_TooMuchData,
               "decode n greater than 256") &&
           ExpectResult(
               leo_ff8xor_decode(
                   64, UINT_MAX, 1, 0,
                   bounded_original, bounded_recovery, bounded_work),
               Leopard_TooMuchData,
               "decode rejects unbounded count before pointer traversal");
}


//------------------------------------------------------------------------------
// End-to-end packed equivalence

struct CodecCase
{
    unsigned OriginalCount;
    unsigned RecoveryCount;
    uint64_t BufferBytes;
    uint32_t Seed;
    const char* Label;
};

struct EncodedData
{
    CodecCase Parameters;
    Buffers PackedOriginal;
    Buffers PackedRecovery;
    Buffers PlaneOriginal;
    Buffers PlaneRecovery;
};

static bool BuildEncodedData(const CodecCase& parameters, EncodedData& encoded)
{
    encoded.Parameters = parameters;
    encoded.PackedOriginal = AllocateBuffers(
        parameters.OriginalCount, parameters.BufferBytes);
    FillRandom(encoded.PackedOriginal, parameters.Seed);

    const unsigned packed_work_count = leo_encode_work_count(
        parameters.OriginalCount, parameters.RecoveryCount);
    Buffers packed_work = AllocateBuffers(
        packed_work_count, parameters.BufferBytes);
    std::vector<const void*> packed_original_ptrs =
        GetConstPointers(encoded.PackedOriginal);
    std::vector<void*> packed_work_ptrs = GetMutablePointers(packed_work);
    const LeopardResult packed_result = leo_encode(
        parameters.BufferBytes,
        parameters.OriginalCount,
        parameters.RecoveryCount,
        packed_work_count,
        packed_original_ptrs.data(),
        packed_work_ptrs.data());
    if (!ExpectResult(packed_result, Leopard_Success, "packed encode"))
        return false;

    encoded.PackedRecovery.assign(
        packed_work.begin(),
        packed_work.begin() + parameters.RecoveryCount);

    encoded.PlaneOriginal = AllocateBuffers(
        parameters.OriginalCount, parameters.BufferBytes);
    for (unsigned index = 0; index < parameters.OriginalCount; ++index)
    {
        PackedToPlane(
            encoded.PackedOriginal[index], encoded.PlaneOriginal[index]);
    }

    const unsigned plane_work_count = leo_ff8xor_encode_work_count(
        parameters.OriginalCount, parameters.RecoveryCount);
    Buffers plane_work = AllocateBuffers(
        plane_work_count, parameters.BufferBytes);
    std::vector<const void*> plane_original_ptrs =
        GetConstPointers(encoded.PlaneOriginal);
    std::vector<void*> plane_work_ptrs = GetMutablePointers(plane_work);
    const LeopardResult plane_result = leo_ff8xor_encode(
        parameters.BufferBytes,
        parameters.OriginalCount,
        parameters.RecoveryCount,
        plane_work_count,
        plane_original_ptrs.data(),
        plane_work_ptrs.data());
    if (!ExpectResult(plane_result, Leopard_Success, "plane encode"))
        return false;

    encoded.PlaneRecovery.assign(
        plane_work.begin(),
        plane_work.begin() + parameters.RecoveryCount);

    for (unsigned index = 0; index < parameters.RecoveryCount; ++index)
    {
        Buffer packed_from_plane;
        PlaneToPacked(encoded.PlaneRecovery[index], packed_from_plane);
        if (!CheckEqual(
                packed_from_plane,
                encoded.PackedRecovery[index],
                parameters.Label,
                index))
            return false;
    }
    return true;
}

static std::vector<unsigned> SelectIndices(
    unsigned total,
    unsigned count,
    uint32_t seed)
{
    std::vector<unsigned> indices(total);
    for (unsigned index = 0; index < total; ++index)
        indices[index] = index;

    for (unsigned remaining = total; remaining > 1; --remaining)
    {
        const unsigned selected = NextRandom(seed) % remaining;
        std::swap(indices[selected], indices[remaining - 1]);
    }
    indices.resize(count);
    std::sort(indices.begin(), indices.end());
    return indices;
}

static bool RunDecodePattern(
    const EncodedData& encoded,
    const std::vector<unsigned>& original_losses,
    const std::vector<unsigned>& recovery_losses,
    LeopardResult expected_result,
    const char* label)
{
    const CodecCase& parameters = encoded.Parameters;
    std::vector<uint8_t> original_missing(parameters.OriginalCount, 0);
    std::vector<uint8_t> recovery_missing(parameters.RecoveryCount, 0);
    for (size_t index = 0; index < original_losses.size(); ++index)
        original_missing[original_losses[index]] = 1;
    for (size_t index = 0; index < recovery_losses.size(); ++index)
        recovery_missing[recovery_losses[index]] = 1;

    std::vector<const void*> packed_original_ptrs =
        GetConstPointers(encoded.PackedOriginal);
    std::vector<const void*> plane_original_ptrs =
        GetConstPointers(encoded.PlaneOriginal);
    std::vector<const void*> packed_recovery_ptrs =
        GetConstPointers(encoded.PackedRecovery);
    std::vector<const void*> plane_recovery_ptrs =
        GetConstPointers(encoded.PlaneRecovery);

    for (unsigned index = 0; index < parameters.OriginalCount; ++index)
    {
        if (original_missing[index])
        {
            packed_original_ptrs[index] = NULL;
            plane_original_ptrs[index] = NULL;
        }
    }
    for (unsigned index = 0; index < parameters.RecoveryCount; ++index)
    {
        if (recovery_missing[index])
        {
            packed_recovery_ptrs[index] = NULL;
            plane_recovery_ptrs[index] = NULL;
        }
    }

    const unsigned packed_work_count = leo_decode_work_count(
        parameters.OriginalCount, parameters.RecoveryCount);
    const unsigned plane_work_count = leo_ff8xor_decode_work_count(
        parameters.OriginalCount, parameters.RecoveryCount);
    Buffers packed_work = AllocateBuffers(
        packed_work_count, parameters.BufferBytes);
    Buffers plane_work = AllocateBuffers(
        plane_work_count, parameters.BufferBytes);
    std::vector<void*> packed_work_ptrs = GetMutablePointers(packed_work);
    std::vector<void*> plane_work_ptrs = GetMutablePointers(plane_work);

    const LeopardResult packed_result = leo_decode(
        parameters.BufferBytes,
        parameters.OriginalCount,
        parameters.RecoveryCount,
        packed_work_count,
        packed_original_ptrs.data(),
        packed_recovery_ptrs.data(),
        packed_work_ptrs.data());
    const LeopardResult plane_result = leo_ff8xor_decode(
        parameters.BufferBytes,
        parameters.OriginalCount,
        parameters.RecoveryCount,
        plane_work_count,
        plane_original_ptrs.data(),
        plane_recovery_ptrs.data(),
        plane_work_ptrs.data());

    if (!ExpectResult(packed_result, expected_result, label) ||
        !ExpectResult(plane_result, expected_result, label))
        return false;
    if (expected_result != Leopard_Success)
        return true;

    for (unsigned index = 0; index < parameters.OriginalCount; ++index)
    {
        if (!original_losses.empty() && !original_missing[index])
            continue;

        if (!CheckEqual(
                packed_work[index],
                encoded.PackedOriginal[index],
                label,
                index) ||
            !CheckEqual(
                plane_work[index],
                encoded.PlaneOriginal[index],
                label,
                index))
            return false;

        Buffer packed_from_plane;
        PlaneToPacked(plane_work[index], packed_from_plane);
        if (!CheckEqual(
                packed_from_plane,
                encoded.PackedOriginal[index],
                label,
                index))
            return false;
    }
    return true;
}

static bool TestDecodePatterns(const EncodedData& encoded)
{
    const CodecCase& p = encoded.Parameters;
    const std::vector<unsigned> none;

    if (!RunDecodePattern(
            encoded, none, none, Leopard_Success, "decode no loss"))
        return false;

    // Keep one fixed loss in every case in addition to the seeded-random
    // subsets below.  This makes the deterministic edge case explicit.
    std::vector<unsigned> one(1, p.OriginalCount / 2);
    if (!RunDecodePattern(
            encoded, one, none, Leopard_Success, "decode one loss"))
        return false;

    const unsigned multiple_count = p.RecoveryCount < 4
        ? p.RecoveryCount
        : 4;
    std::vector<unsigned> multiple = SelectIndices(
        p.OriginalCount, multiple_count, p.Seed ^ 0x22222222U);
    if (multiple_count > 1 &&
        !RunDecodePattern(
            encoded, multiple, none,
            Leopard_Success, "decode multiple losses"))
        return false;

    if (multiple_count != p.RecoveryCount)
    {
        std::vector<unsigned> maximum = SelectIndices(
            p.OriginalCount,
            p.RecoveryCount,
            p.Seed ^ 0x33333333U);
        if (!RunDecodePattern(
                encoded, maximum, none,
                Leopard_Success, "decode maximum losses"))
            return false;
    }

    if (p.RecoveryCount > 1)
    {
        const unsigned mixed_original_count = p.RecoveryCount - 1 < 3
            ? p.RecoveryCount - 1
            : 3;
        const unsigned mixed_recovery_count = p.RecoveryCount -
            mixed_original_count < 2
            ? p.RecoveryCount - mixed_original_count
            : 2;
        std::vector<unsigned> mixed_original = SelectIndices(
            p.OriginalCount,
            mixed_original_count,
            p.Seed ^ 0x44444444U);
        std::vector<unsigned> mixed_recovery = SelectIndices(
            p.RecoveryCount,
            mixed_recovery_count,
            p.Seed ^ 0x55555555U);
        if (!RunDecodePattern(
                encoded,
                mixed_original,
                mixed_recovery,
                Leopard_Success,
                "decode mixed original/recovery losses"))
            return false;
    }

    const unsigned insufficient_original_count = p.RecoveryCount < 2
        ? 1
        : 2;
    const unsigned recovery_to_keep = insufficient_original_count - 1;
    std::vector<unsigned> insufficient_original = SelectIndices(
        p.OriginalCount,
        insufficient_original_count,
        p.Seed ^ 0x66666666U);
    std::vector<unsigned> insufficient_recovery = SelectIndices(
        p.RecoveryCount,
        p.RecoveryCount - recovery_to_keep,
        p.Seed ^ 0x77777777U);
    return RunDecodePattern(
        encoded,
        insufficient_original,
        insufficient_recovery,
        Leopard_NeedMoreData,
        "decode insufficient recovery");
}

static bool TestPackedEquivalence()
{
    static const CodecCase kCases[] = {
        {   1,   1,    64, 0x10010001U, "encode k=1 r=1 bytes=64" },
        {   8,   1,  1024, 0x10080001U, "encode k=8 r=1 bytes=1024" },
        {   4,   2,    64, 0x20040002U, "encode k=4 r=2 bytes=64" },
        {   7,   3,  1024, 0x30070003U, "encode k=7 r=3 bytes=1024" },
        {   8,   2,  1024, 0x20080002U, "encode k=8 r=2 bytes=1024" },
        {  16,   4,  4096, 0x40160004U, "encode k=16 r=4 bytes=4096" },
        {  17,   5,  4096, 0x50170005U, "encode k=17 r=5 bytes=4096" },
        {  32,   8, 65536, 0x80320008U, "encode k=32 r=8 bytes=65536" },
        {  64,  16,  1024, 0x10640010U, "encode k=64 r=16 bytes=1024" },
        { 128,  32,  4096, 0x20800020U, "encode k=128 r=32 bytes=4096" },
        { 128, 128,  1024, 0x80800080U, "encode k=128 r=128 bytes=1024" }
    };

    for (unsigned index = 0;
         index < sizeof(kCases) / sizeof(kCases[0]);
         ++index)
    {
        EncodedData encoded;
        if (!BuildEncodedData(kCases[index], encoded) ||
            !TestDecodePatterns(encoded))
            return false;
    }
    return true;
}

static bool TestForcedLocatorShifts()
{
    const CodecCase parameters = {
        8, 4, 64, 0x25510ca7U, "forced locator shift"
    };
    EncodedData encoded;
    if (!BuildEncodedData(parameters, encoded))
        return false;

    // Two unavailable originals and one unavailable recovery exercise all
    // three locator-scaling sites: available inputs, unavailable inputs, and
    // final inverse scaling of recovered originals.
    static const unsigned kOriginalLosses[] = { 1, 6 };
    static const unsigned kRecoveryLosses[] = { 2 };
    std::vector<const void*> original_ptrs =
        GetConstPointers(encoded.PlaneOriginal);
    std::vector<const void*> recovery_ptrs =
        GetConstPointers(encoded.PlaneRecovery);
    for (unsigned index = 0;
         index < sizeof(kOriginalLosses) / sizeof(kOriginalLosses[0]);
         ++index)
        original_ptrs[kOriginalLosses[index]] = NULL;
    for (unsigned index = 0;
         index < sizeof(kRecoveryLosses) / sizeof(kRecoveryLosses[0]);
         ++index)
        recovery_ptrs[kRecoveryLosses[index]] = NULL;

    const unsigned work_count = leo_ff8xor_decode_work_count(
        parameters.OriginalCount, parameters.RecoveryCount);
    Buffers work = AllocateBuffers(work_count, parameters.BufferBytes);
    std::vector<void*> work_ptrs = GetMutablePointers(work);

    bool success = true;
    for (unsigned shift = 0; shift < kFieldModulus; ++shift)
    {
        for (unsigned index = 0; index < work_count; ++index)
        {
            std::fill(
                work[index].begin(),
                work[index].end(),
                static_cast<uint8_t>(0x5aU ^ shift));
        }

        leopard::ff8xor::SetLocatorShiftForTesting(
            static_cast<int>(shift));
        const LeopardResult result = leo_ff8xor_decode(
            parameters.BufferBytes,
            parameters.OriginalCount,
            parameters.RecoveryCount,
            work_count,
            original_ptrs.data(),
            recovery_ptrs.data(),
            work_ptrs.data());
        if (!ExpectResult(result, Leopard_Success, parameters.Label))
        {
            success = false;
            break;
        }
        if (leopard::ff8xor::GetLastLocatorShiftForTesting() != shift)
        {
            fprintf(stderr,
                "forced locator shift mismatch: requested=%u observed=%u\n",
                shift,
                leopard::ff8xor::GetLastLocatorShiftForTesting());
            success = false;
            break;
        }

        for (unsigned loss_index = 0;
             loss_index <
                sizeof(kOriginalLosses) / sizeof(kOriginalLosses[0]);
             ++loss_index)
        {
            const unsigned original_index = kOriginalLosses[loss_index];
            if (!CheckEqual(
                    work[original_index],
                    encoded.PlaneOriginal[original_index],
                    parameters.Label,
                    original_index))
            {
                fprintf(stderr, "failed forced locator shift %u\n", shift);
                success = false;
                break;
            }

            Buffer packed_from_plane;
            PlaneToPacked(work[original_index], packed_from_plane);
            if (!CheckEqual(
                    packed_from_plane,
                    encoded.PackedOriginal[original_index],
                    parameters.Label,
                    original_index))
            {
                fprintf(stderr, "failed forced locator shift %u\n", shift);
                success = false;
                break;
            }
        }
        if (!success)
            break;
    }

    leopard::ff8xor::SetLocatorShiftForTesting(-1);
    return success;
}


//------------------------------------------------------------------------------
// Native plane-layout round trips

static bool RunNativeRoundTrip(const CodecCase& parameters)
{
    Buffers original = AllocateBuffers(
        parameters.OriginalCount, parameters.BufferBytes);
    FillRandom(original, parameters.Seed);
    const unsigned encode_work_count = leo_ff8xor_encode_work_count(
        parameters.OriginalCount, parameters.RecoveryCount);
    Buffers encode_work = AllocateBuffers(
        encode_work_count, parameters.BufferBytes);
    std::vector<const void*> original_ptrs = GetConstPointers(original);
    std::vector<void*> encode_work_ptrs = GetMutablePointers(encode_work);
    if (!ExpectResult(
            leo_ff8xor_encode(
                parameters.BufferBytes,
                parameters.OriginalCount,
                parameters.RecoveryCount,
                encode_work_count,
                original_ptrs.data(),
                encode_work_ptrs.data()),
            Leopard_Success,
            "native encode"))
        return false;

    Buffers recovery(
        encode_work.begin(),
        encode_work.begin() + parameters.RecoveryCount);
    std::vector<const void*> recovery_ptrs = GetConstPointers(recovery);

    const unsigned loss_count = parameters.RecoveryCount < 4
        ? parameters.RecoveryCount
        : 4;
    std::vector<unsigned> losses = SelectIndices(
        parameters.OriginalCount,
        loss_count,
        parameters.Seed ^ 0xa5a5a5a5U);
    for (size_t index = 0; index < losses.size(); ++index)
        original_ptrs[losses[index]] = NULL;

    if (parameters.RecoveryCount > loss_count)
    {
        const unsigned recovery_losses = parameters.RecoveryCount - loss_count < 2
            ? parameters.RecoveryCount - loss_count
            : 2;
        std::vector<unsigned> missing = SelectIndices(
            parameters.RecoveryCount,
            recovery_losses,
            parameters.Seed ^ 0x5a5a5a5aU);
        for (size_t index = 0; index < missing.size(); ++index)
            recovery_ptrs[missing[index]] = NULL;
    }

    const unsigned decode_work_count = leo_ff8xor_decode_work_count(
        parameters.OriginalCount, parameters.RecoveryCount);
    Buffers decode_work = AllocateBuffers(
        decode_work_count, parameters.BufferBytes);
    std::vector<void*> decode_work_ptrs = GetMutablePointers(decode_work);
    if (!ExpectResult(
            leo_ff8xor_decode(
                parameters.BufferBytes,
                parameters.OriginalCount,
                parameters.RecoveryCount,
                decode_work_count,
                original_ptrs.data(),
                recovery_ptrs.data(),
                decode_work_ptrs.data()),
            Leopard_Success,
            parameters.Label))
        return false;

    for (size_t index = 0; index < losses.size(); ++index)
    {
        const unsigned original_index = losses[index];
        if (!CheckEqual(
                decode_work[original_index],
                original[original_index],
                parameters.Label,
                original_index))
            return false;
    }
    return true;
}

static bool TestNativeRoundTrips()
{
    static const CodecCase kCases[] = {
        {  4, 2,   64, 0x64040002U, "native k=4 r=2 bytes=64" },
        {  7, 3,  192, 0x19070003U, "native k=7 r=3 bytes=192" },
        {  8, 2,  192, 0x19080002U, "native k=8 r=2 bytes=192" },
        { 16, 4,  320, 0x32160004U, "native k=16 r=4 bytes=320" },
        { 32, 8, 1024, 0x10320008U, "native k=32 r=8 bytes=1024" }
    };
    for (unsigned index = 0;
         index < sizeof(kCases) / sizeof(kCases[0]);
         ++index)
    {
        if (!RunNativeRoundTrip(kCases[index]))
            return false;
    }
    return true;
}


//------------------------------------------------------------------------------
// Direct portable/SIMD kernel checks

static bool CheckMultiplyKernel(uint64_t buffer_bytes, unsigned coefficient)
{
    uint32_t state = static_cast<uint32_t>(
        0xc001d00dU ^ coefficient ^ static_cast<unsigned>(buffer_bytes));
    Buffer packed_source(static_cast<size_t>(buffer_bytes));
    Buffer packed_expected(static_cast<size_t>(buffer_bytes));
    FillRandom(packed_source, state);
    for (size_t index = 0; index < packed_source.size(); ++index)
    {
        packed_expected[index] = ReferenceMultiply(
            packed_source[index], static_cast<uint8_t>(coefficient));
    }

    Buffer plane_source;
    Buffer plane_expected;
    PackedToPlane(packed_source, plane_source);
    PackedToPlane(packed_expected, plane_expected);

    Buffer out_of_place(static_cast<size_t>(buffer_bytes), 0xa5);
    leopard::ff8xor::MultiplyBuffer(
        buffer_bytes,
        out_of_place.data(),
        plane_source.data(),
        static_cast<uint8_t>(coefficient));
    if (!CheckEqual(
            out_of_place, plane_expected, "direct multiply", coefficient))
        return false;

    Buffer in_place = plane_source;
    leopard::ff8xor::MultiplyBuffer(
        buffer_bytes,
        in_place.data(),
        in_place.data(),
        static_cast<uint8_t>(coefficient));
    return CheckEqual(
        in_place, plane_expected, "in-place direct multiply", coefficient);
}

static bool CheckIdentityMultiplyFastPath(
    uint64_t buffer_bytes,
    unsigned coefficient)
{
    const size_t bytes = static_cast<size_t>(buffer_bytes);
    const size_t overlap = 13;
    uint32_t state = static_cast<uint32_t>(
        0x1de1717eU ^ coefficient ^ static_cast<unsigned>(buffer_bytes));

    // A nonzero byte count with null in-place pointers is deliberate:  The
    // identity path is required to return before attempting any payload load
    // or store.  This also catches accidental dispatch to the generated
    // identity circuit, whose zero gates still surround payload I/O.
    leopard::ff8xor::MultiplyBuffer(
        buffer_bytes,
        static_cast<void*>(NULL),
        static_cast<const void*>(NULL),
        static_cast<uint8_t>(coefficient));

    Buffer forward(bytes + overlap);
    FillRandom(forward, state);
    Buffer expected_source(forward.begin(), forward.begin() + bytes);
    leopard::ff8xor::MultiplyBuffer(
        buffer_bytes,
        forward.data() + overlap,
        forward.data(),
        static_cast<uint8_t>(coefficient));
    if (!std::equal(
            expected_source.begin(),
            expected_source.end(),
            forward.begin() + overlap))
    {
        fprintf(stderr,
            "overlapping forward identity copy mismatch: bytes=%llu log=%u\n",
            static_cast<unsigned long long>(buffer_bytes), coefficient);
        return false;
    }

    Buffer backward(bytes + overlap);
    FillRandom(backward, state);
    expected_source.assign(backward.begin() + overlap, backward.end());
    leopard::ff8xor::MultiplyBuffer(
        buffer_bytes,
        backward.data(),
        backward.data() + overlap,
        static_cast<uint8_t>(coefficient));
    if (!std::equal(
            expected_source.begin(),
            expected_source.end(),
            backward.begin()))
    {
        fprintf(stderr,
            "overlapping backward identity copy mismatch: bytes=%llu log=%u\n",
            static_cast<unsigned long long>(buffer_bytes), coefficient);
        return false;
    }

    return true;
}

static bool CheckButterflyKernel(uint64_t buffer_bytes, unsigned skew)
{
    uint32_t state = static_cast<uint32_t>(
        0xb17ef17eU ^ skew ^ static_cast<unsigned>(buffer_bytes));
    Buffer packed_x(static_cast<size_t>(buffer_bytes));
    Buffer packed_y(static_cast<size_t>(buffer_bytes));
    Buffer expected_fft_x(static_cast<size_t>(buffer_bytes));
    Buffer expected_fft_y(static_cast<size_t>(buffer_bytes));
    Buffer expected_ifft_x(static_cast<size_t>(buffer_bytes));
    Buffer expected_ifft_y(static_cast<size_t>(buffer_bytes));
    FillRandom(packed_x, state);
    FillRandom(packed_y, state);

    for (size_t index = 0; index < packed_x.size(); ++index)
    {
        const uint16_t input = static_cast<uint16_t>(
            packed_x[index] | (static_cast<uint16_t>(packed_y[index]) << 8));
        const uint16_t fft = ReferenceFFT(input, skew);
        const uint16_t ifft = ReferenceIFFT(input, skew);
        expected_fft_x[index] = static_cast<uint8_t>(fft);
        expected_fft_y[index] = static_cast<uint8_t>(fft >> 8);
        expected_ifft_x[index] = static_cast<uint8_t>(ifft);
        expected_ifft_y[index] = static_cast<uint8_t>(ifft >> 8);
    }

    Buffer plane_x;
    Buffer plane_y;
    Buffer expected_plane_x;
    Buffer expected_plane_y;
    PackedToPlane(packed_x, plane_x);
    PackedToPlane(packed_y, plane_y);
    PackedToPlane(expected_fft_x, expected_plane_x);
    PackedToPlane(expected_fft_y, expected_plane_y);
    leopard::ff8xor::FFTButterflyBuffer(
        buffer_bytes,
        plane_x.data(),
        plane_y.data(),
        static_cast<uint8_t>(skew));
    if (!CheckEqual(plane_x, expected_plane_x, "direct FFT x", skew) ||
        !CheckEqual(plane_y, expected_plane_y, "direct FFT y", skew))
        return false;

    leopard::ff8xor::IFFTButterflyBuffer(
        buffer_bytes,
        plane_x.data(),
        plane_y.data(),
        static_cast<uint8_t>(skew));
    Buffer original_plane_x;
    Buffer original_plane_y;
    PackedToPlane(packed_x, original_plane_x);
    PackedToPlane(packed_y, original_plane_y);
    if (!CheckEqual(plane_x, original_plane_x, "FFT/IFFT inverse x", skew) ||
        !CheckEqual(plane_y, original_plane_y, "FFT/IFFT inverse y", skew))
        return false;

    PackedToPlane(packed_x, plane_x);
    PackedToPlane(packed_y, plane_y);
    PackedToPlane(expected_ifft_x, expected_plane_x);
    PackedToPlane(expected_ifft_y, expected_plane_y);
    leopard::ff8xor::IFFTButterflyBuffer(
        buffer_bytes,
        plane_x.data(),
        plane_y.data(),
        static_cast<uint8_t>(skew));
    return CheckEqual(plane_x, expected_plane_x, "direct IFFT x", skew) &&
           CheckEqual(plane_y, expected_plane_y, "direct IFFT y", skew);
}

static bool TestKernelModes()
{
    typedef leopard::ff8xor::KernelMode KernelMode;
    static const KernelMode kModes[] = {
        KernelMode::Portable,
        KernelMode::Simd128,
        KernelMode::Avx2,
        KernelMode::Avx512VL,
        KernelMode::Avx512Zmm
    };
    static const char* kModeNames[] = {
        "portable uint64 XOR circuits",
        "128-bit SIMD XOR circuits",
        "AVX2 XOR circuits",
        "AVX-512VL YMM XOR circuits",
        "AVX-512 ZMM XOR circuits"
    };
    // The larger sizes force both AVX-512 widths, portable and 128-bit tails,
    // and an AVX2 tail after a 512-bit ZMM chunk.
    static const uint64_t kSizes[] = {
        64, 128, 192, 256, 320, 512, 576, 640, 768, 1024, 1536
    };
    static const unsigned kCoefficients[] = { 0, 1, 51, 254, 255 };
    static const unsigned kSkews[] = { 0, 51, 254, 255 };
    const KernelMode saved_mode = leopard::ff8xor::GetKernelMode();

    for (unsigned mode_index = 0;
         mode_index < sizeof(kModes) / sizeof(kModes[0]);
         ++mode_index)
    {
        const KernelMode mode = kModes[mode_index];
        if (!leopard::ff8xor::IsKernelModeAvailable(mode))
        {
            // An unavailable forced mode must resolve to a safe baseline mode,
            // never remain selected and risk an illegal instruction.
            leopard::ff8xor::SetKernelMode(mode);
            if (leopard::ff8xor::GetActiveKernelMode() == mode ||
                !CheckMultiplyKernel(512, 51) ||
                !CheckButterflyKernel(512, 51))
            {
                fprintf(stderr,
                    "unavailable forced kernel mode did not fall back safely\n");
                leopard::ff8xor::SetKernelMode(saved_mode);
                return false;
            }
            continue;
        }
        leopard::ff8xor::SetKernelMode(mode);
        if (leopard::ff8xor::GetActiveKernelMode() != mode)
        {
            fprintf(stderr, "requested kernel mode did not become active\n");
            leopard::ff8xor::SetKernelMode(saved_mode);
            return false;
        }
        if (strcmp(
                leopard::ff8xor::GetKernelBackendName(),
                kModeNames[mode_index]) != 0)
        {
            fprintf(stderr, "active kernel mode name mismatch\n");
            leopard::ff8xor::SetKernelMode(saved_mode);
            return false;
        }

        for (unsigned size_index = 0;
             size_index < sizeof(kSizes) / sizeof(kSizes[0]);
             ++size_index)
        {
            for (unsigned coefficient_index = 0;
                 coefficient_index <
                    sizeof(kCoefficients) / sizeof(kCoefficients[0]);
                 ++coefficient_index)
            {
                if (!CheckMultiplyKernel(
                        kSizes[size_index],
                        kCoefficients[coefficient_index]))
                {
                    leopard::ff8xor::SetKernelMode(saved_mode);
                    return false;
                }

                const unsigned coefficient =
                    kCoefficients[coefficient_index];
                if ((coefficient == 0 || coefficient == kFieldModulus) &&
                    !CheckIdentityMultiplyFastPath(
                        kSizes[size_index], coefficient))
                {
                    leopard::ff8xor::SetKernelMode(saved_mode);
                    return false;
                }
            }

            for (unsigned skew_index = 0;
                 skew_index < sizeof(kSkews) / sizeof(kSkews[0]);
                 ++skew_index)
            {
                if (!CheckButterflyKernel(
                        kSizes[size_index],
                        kSkews[skew_index]))
                {
                    leopard::ff8xor::SetKernelMode(saved_mode);
                    return false;
                }
            }
        }
    }

    leopard::ff8xor::SetKernelMode(saved_mode);
    return true;
}

static bool TestAllAVX512Specializations()
{
    typedef leopard::ff8xor::KernelMode KernelMode;
    static const KernelMode kModes[] = {
        KernelMode::Avx512VL,
        KernelMode::Avx512Zmm
    };
    const KernelMode saved_mode = leopard::ff8xor::GetKernelMode();

    for (unsigned mode_index = 0;
         mode_index < sizeof(kModes) / sizeof(kModes[0]);
         ++mode_index)
    {
        const KernelMode mode = kModes[mode_index];
        if (!leopard::ff8xor::IsKernelModeAvailable(mode))
            continue;

        leopard::ff8xor::SetKernelMode(mode);

        // 1536 bytes gives 192 bytes per plane: Six YMM chunks or three ZMM
        // chunks.  The size sweep above separately exercises AVX2, 128-bit,
        // and uint64_t tails.  Identity multipliers 0 and 255 and butterfly
        // sentinel 255 are intentionally handled by the separately tested
        // no-op/contiguous-XOR fast paths.
        for (unsigned coefficient = 1;
             coefficient < kFieldModulus;
             ++coefficient)
        {
            if (!CheckMultiplyKernel(1536, coefficient))
            {
                fprintf(stderr,
                    "AVX-512 multiply specialization failed: mode=%u log=%u\n",
                    static_cast<unsigned>(mode), coefficient);
                leopard::ff8xor::SetKernelMode(saved_mode);
                return false;
            }
        }

        for (unsigned skew = 0; skew < kFieldModulus; ++skew)
        {
            if (!CheckButterflyKernel(1536, skew))
            {
                fprintf(stderr,
                    "AVX-512 butterfly specialization failed: mode=%u skew=%u\n",
                    static_cast<unsigned>(mode), skew);
                leopard::ff8xor::SetKernelMode(saved_mode);
                return false;
            }
        }

        // Exercise the complete encode/decode control flow with the forced
        // mode rather than validating only isolated whole-buffer operations.
        CodecCase codec_case = {
            16, 4, 1536,
            static_cast<uint32_t>(0xa5120000U + mode_index),
            "forced AVX-512 packed equivalence"
        };
        EncodedData encoded;
        if (!BuildEncodedData(codec_case, encoded) ||
            !TestDecodePatterns(encoded))
        {
            leopard::ff8xor::SetKernelMode(saved_mode);
            return false;
        }
    }

    leopard::ff8xor::SetKernelMode(saved_mode);
    return true;
}

static bool TestAVX2FeaturePredicate()
{
    const uint32_t osxsave = UINT32_C(1) << 27;
    const uint32_t avx = UINT32_C(1) << 28;
    const uint32_t avx2 = UINT32_C(1) << 5;
    const uint32_t leaf1 = osxsave | avx;

    if (!leopard::IsAVX2Supported(7, leaf1, avx2, UINT64_C(0x6)) ||
        !leopard::IsAVX2Supported(0xffffffffU, 0xffffffffU,
            0xffffffffU, UINT64_MAX) ||
        leopard::IsAVX2Supported(6, leaf1, avx2, UINT64_C(0x6)) ||
        leopard::IsAVX2Supported(7, avx, avx2, UINT64_C(0x6)) ||
        leopard::IsAVX2Supported(7, osxsave, avx2, UINT64_C(0x6)) ||
        leopard::IsAVX2Supported(7, leaf1, 0, UINT64_C(0x6)) ||
        leopard::IsAVX2Supported(7, leaf1, avx2, UINT64_C(0x0)) ||
        leopard::IsAVX2Supported(7, leaf1, avx2, UINT64_C(0x2)) ||
        leopard::IsAVX2Supported(7, leaf1, avx2, UINT64_C(0x4)))
    {
        fprintf(stderr, "AVX2 CPU/OS feature predicate failed\n");
        return false;
    }
    return true;
}

static bool TestNextPow2Boundaries()
{
    if (leopard::NextPow2(0) != 0 ||
        leopard::NextPow2(1) != 1 ||
        leopard::NextPow2(2) != 2 ||
        leopard::NextPow2(3) != 4 ||
        leopard::NextPow2(UINT32_C(0x80000000)) !=
            UINT32_C(0x80000000) ||
        leopard::NextPow2(UINT32_C(0x80000001)) != 0 ||
        leopard::NextPow2(UINT32_MAX) != 0)
    {
        fprintf(stderr, "NextPow2 boundary handling failed\n");
        return false;
    }
    return true;
}

} // namespace

int main()
{
    InitializeReferenceField();

    if (!TestMultiplicationCircuits() ||
        !TestButterflyCircuits() ||
        !TestTransposeHelpers() ||
        !TestWorkCounts() ||
        !TestCallInitialize() ||
        !TestAVX2FeaturePredicate() ||
        !TestNextPow2Boundaries())
        return 1;

    if (!ExpectResult(
            static_cast<LeopardResult>(leo_init()),
        Leopard_Success,
        "leo_init") ||
        !TestAPIValidation() ||
        !TestKernelModes() ||
        !TestAllAVX512Specializations() ||
        !TestPackedEquivalence() ||
        !TestForcedLocatorShifts() ||
        !TestNativeRoundTrips())
        return 1;

    printf("FF8 XOR circuit tests passed\n");
    printf("avx512vl=%s avx512zmm=%s\n",
        leopard::ff8xor::IsKernelModeAvailable(
            leopard::ff8xor::KernelMode::Avx512VL) ? "available" : "unavailable",
        leopard::ff8xor::IsKernelModeAvailable(
            leopard::ff8xor::KernelMode::Avx512Zmm) ? "available" : "unavailable");
    printf("checksum=%s\n", leopard::ff8xor::generated::kCircuitChecksum);
    printf("multiply_gates=%u..%u avg=%.6f\n",
        leopard::ff8xor::generated::kMultiplyMinGateCount,
        leopard::ff8xor::generated::kMultiplyMaxGateCount,
        leopard::ff8xor::generated::kMultiplyAverageGateCount);
    printf("fft_gates=%u..%u avg=%.6f ifft_gates=%u..%u avg=%.6f\n",
        leopard::ff8xor::generated::kFFTMinGateCount,
        leopard::ff8xor::generated::kFFTMaxGateCount,
        leopard::ff8xor::generated::kFFTAverageGateCount,
        leopard::ff8xor::generated::kIFFTMinGateCount,
        leopard::ff8xor::generated::kIFFTMaxGateCount,
        leopard::ff8xor::generated::kIFFTAverageGateCount);
    return 0;
}
