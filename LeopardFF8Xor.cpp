/*
    Copyright (c) 2017 Christopher A. Taylor.  All rights reserved.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are met:

    * Redistributions of source code must retain the above copyright notice,
      this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright notice,
      this list of conditions and the following disclaimer in the documentation
      and/or other materials provided with the distribution.
    * Neither the name of Leopard-RS nor the names of its contributors may be
      used to endorse or promote products derived from this software without
      specific prior written permission.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
    AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
    IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
    ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
    LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
    CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
    SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
    INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
    CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
    ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
    POSSIBILITY OF SUCH DAMAGE.
*/

#include "LeopardFF8Xor.h"

#ifdef LEO_HAS_FF8

#include <atomic>
#include <limits.h>
#include <string.h>

#if !defined(LEO_TARGET_MOBILE) || defined(LEO_TRY_NEON) || defined(LEO_USE_SSE2NEON)
    #define LEO_FF8XOR_HAS_SIMD128
#endif

namespace leopard { namespace ff8xor {


// The generated circuit expressions use this overload set.  Keeping the
// result in the destination operand lets x86 compilers select three-operand
// VPXOR without introducing an additional vector temporary.
static LEO_FORCE_INLINE uint64_t XorValue(uint64_t x, uint64_t y)
{
    return x ^ y;
}

#ifdef LEO_FF8XOR_HAS_SIMD128
static LEO_FORCE_INLINE LEO_M128 XorValue(LEO_M128 x, LEO_M128 y)
{
#if defined(LEO_TRY_NEON) && !defined(LEO_USE_SSE2NEON)
    return veorq_u8(x, y);
#else
    return _mm_xor_si128(x, y);
#endif
}
#endif

#ifdef LEO_TRY_AVX2
static LEO_FORCE_INLINE LEO_M256 XorValue(LEO_M256 x, LEO_M256 y)
{
    return _mm256_xor_si256(x, y);
}
#endif


}} // namespace leopard::ff8xor

#include "generated/LeopardFF8XorCircuits.inl"

namespace leopard { namespace ff8xor {


//------------------------------------------------------------------------------
// Plane value loads and stores

template <typename Value>
struct ValueIO;

struct PortableTag {};

template <>
struct ValueIO<PortableTag>
{
    typedef uint64_t Value;
    static const unsigned kBytes = 8;

    static LEO_FORCE_INLINE Value Load(const uint8_t* source)
    {
        Value value;
        memcpy(&value, source, sizeof(value));
        return value;
    }

    static LEO_FORCE_INLINE void Store(uint8_t* destination, Value value)
    {
        memcpy(destination, &value, sizeof(value));
    }
};

#ifdef LEO_FF8XOR_HAS_SIMD128
struct Simd128Tag {};

template <>
struct ValueIO<Simd128Tag>
{
    typedef LEO_M128 Value;
    static const unsigned kBytes = 16;

    static LEO_FORCE_INLINE Value Load(const uint8_t* source)
    {
#if defined(LEO_TRY_NEON) && !defined(LEO_USE_SSE2NEON)
        return vld1q_u8(source);
#else
        return _mm_loadu_si128(reinterpret_cast<const LEO_M128*>(source));
#endif
    }

    static LEO_FORCE_INLINE void Store(uint8_t* destination, Value value)
    {
#if defined(LEO_TRY_NEON) && !defined(LEO_USE_SSE2NEON)
        vst1q_u8(destination, value);
#else
        _mm_storeu_si128(reinterpret_cast<LEO_M128*>(destination), value);
#endif
    }
};
#endif

#ifdef LEO_TRY_AVX2
struct Avx2Tag {};

template <>
struct ValueIO<Avx2Tag>
{
    typedef LEO_M256 Value;
    static const unsigned kBytes = 32;

    static LEO_FORCE_INLINE Value Load(const uint8_t* source)
    {
        return _mm256_loadu_si256(reinterpret_cast<const LEO_M256*>(source));
    }

    static LEO_FORCE_INLINE void Store(uint8_t* destination, Value value)
    {
        _mm256_storeu_si256(reinterpret_cast<LEO_M256*>(destination), value);
    }
};
#endif


//------------------------------------------------------------------------------
// Kernel selection

// Test inspection and overrides must not introduce shared mutable state into
// otherwise independent codec calls.
static std::atomic<KernelMode> SelectedKernelMode(KernelMode::Auto);
static thread_local unsigned LastLocatorShift = 0;
static thread_local int LocatorShiftOverride = -1;

bool IsKernelModeAvailable(KernelMode mode)
{
    switch (mode)
    {
    case KernelMode::Auto:
    case KernelMode::Portable:
        return true;

    case KernelMode::Simd128:
#ifdef LEO_FF8XOR_HAS_SIMD128
        return true;
#else
        return false;
#endif

    case KernelMode::Avx2:
#ifdef LEO_TRY_AVX2
        return CpuHasAVX2;
#else
        return false;
#endif
    }
    return false;
}

static KernelMode ResolveKernelMode()
{
    const KernelMode selected =
        SelectedKernelMode.load(std::memory_order_relaxed);

    if (selected == KernelMode::Portable)
        return KernelMode::Portable;

    if (selected == KernelMode::Avx2)
    {
#ifdef LEO_TRY_AVX2
        if (CpuHasAVX2)
            return KernelMode::Avx2;
#endif
#ifdef LEO_FF8XOR_HAS_SIMD128
        return KernelMode::Simd128;
#else
        return KernelMode::Portable;
#endif
    }

    if (selected == KernelMode::Simd128)
    {
#ifdef LEO_FF8XOR_HAS_SIMD128
        return KernelMode::Simd128;
#else
        return KernelMode::Portable;
#endif
    }

#ifdef LEO_TRY_AVX2
    if (CpuHasAVX2)
        return KernelMode::Avx2;
#endif
#ifdef LEO_FF8XOR_HAS_SIMD128
    return KernelMode::Simd128;
#else
    return KernelMode::Portable;
#endif
}

void SetKernelMode(KernelMode mode)
{
    SelectedKernelMode.store(mode, std::memory_order_relaxed);
}

KernelMode GetKernelMode()
{
    return SelectedKernelMode.load(std::memory_order_relaxed);
}

KernelMode GetActiveKernelMode()
{
    return ResolveKernelMode();
}

const char* GetKernelBackendName()
{
    switch (ResolveKernelMode())
    {
    case KernelMode::Avx2: return "AVX2 XOR circuits";
    case KernelMode::Simd128: return "128-bit SIMD XOR circuits";
    case KernelMode::Portable: return "portable uint64 XOR circuits";
    case KernelMode::Auto: break;
    }
    return "unknown XOR circuit backend";
}

unsigned GetLastLocatorShiftForTesting()
{
    return LastLocatorShift;
}

void SetLocatorShiftForTesting(int shift)
{
    LocatorShiftOverride = shift >= 0 && shift < kModulus ? shift : -1;
}


//------------------------------------------------------------------------------
// Named-register circuit chunks

template <unsigned Coefficient, typename ValueTag>
static LEO_FORCE_INLINE void MultiplyChunk(
    uint8_t* destination,
    const uint8_t* source,
    uint64_t plane_bytes,
    uint64_t offset)
{
    typedef typename ValueIO<ValueTag>::Value Value;
    Value x0 = ValueIO<ValueTag>::Load(source + offset);
    Value x1 = ValueIO<ValueTag>::Load(source + plane_bytes + offset);
    Value x2 = ValueIO<ValueTag>::Load(source + plane_bytes * 2 + offset);
    Value x3 = ValueIO<ValueTag>::Load(source + plane_bytes * 3 + offset);
    Value x4 = ValueIO<ValueTag>::Load(source + plane_bytes * 4 + offset);
    Value x5 = ValueIO<ValueTag>::Load(source + plane_bytes * 5 + offset);
    Value x6 = ValueIO<ValueTag>::Load(source + plane_bytes * 6 + offset);
    Value x7 = ValueIO<ValueTag>::Load(source + plane_bytes * 7 + offset);

    generated::MultiplyCircuit<Coefficient>::Apply(
        x0, x1, x2, x3, x4, x5, x6, x7);

    ValueIO<ValueTag>::Store(destination + offset, x0);
    ValueIO<ValueTag>::Store(destination + plane_bytes + offset, x1);
    ValueIO<ValueTag>::Store(destination + plane_bytes * 2 + offset, x2);
    ValueIO<ValueTag>::Store(destination + plane_bytes * 3 + offset, x3);
    ValueIO<ValueTag>::Store(destination + plane_bytes * 4 + offset, x4);
    ValueIO<ValueTag>::Store(destination + plane_bytes * 5 + offset, x5);
    ValueIO<ValueTag>::Store(destination + plane_bytes * 6 + offset, x6);
    ValueIO<ValueTag>::Store(destination + plane_bytes * 7 + offset, x7);
}

template <unsigned Skew, typename ValueTag, bool Inverse>
static LEO_FORCE_INLINE void ButterflyChunk(
    uint8_t* x_buffer,
    uint8_t* y_buffer,
    uint64_t plane_bytes,
    uint64_t offset)
{
    typedef typename ValueIO<ValueTag>::Value Value;
    Value x0 = ValueIO<ValueTag>::Load(x_buffer + offset);
    Value x1 = ValueIO<ValueTag>::Load(x_buffer + plane_bytes + offset);
    Value x2 = ValueIO<ValueTag>::Load(x_buffer + plane_bytes * 2 + offset);
    Value x3 = ValueIO<ValueTag>::Load(x_buffer + plane_bytes * 3 + offset);
    Value x4 = ValueIO<ValueTag>::Load(x_buffer + plane_bytes * 4 + offset);
    Value x5 = ValueIO<ValueTag>::Load(x_buffer + plane_bytes * 5 + offset);
    Value x6 = ValueIO<ValueTag>::Load(x_buffer + plane_bytes * 6 + offset);
    Value x7 = ValueIO<ValueTag>::Load(x_buffer + plane_bytes * 7 + offset);
    Value y0 = ValueIO<ValueTag>::Load(y_buffer + offset);
    Value y1 = ValueIO<ValueTag>::Load(y_buffer + plane_bytes + offset);
    Value y2 = ValueIO<ValueTag>::Load(y_buffer + plane_bytes * 2 + offset);
    Value y3 = ValueIO<ValueTag>::Load(y_buffer + plane_bytes * 3 + offset);
    Value y4 = ValueIO<ValueTag>::Load(y_buffer + plane_bytes * 4 + offset);
    Value y5 = ValueIO<ValueTag>::Load(y_buffer + plane_bytes * 5 + offset);
    Value y6 = ValueIO<ValueTag>::Load(y_buffer + plane_bytes * 6 + offset);
    Value y7 = ValueIO<ValueTag>::Load(y_buffer + plane_bytes * 7 + offset);

    if (Inverse)
    {
        generated::IFFTCircuit<Skew>::Apply(
            x0, x1, x2, x3, x4, x5, x6, x7,
            y0, y1, y2, y3, y4, y5, y6, y7);
    }
    else
    {
        generated::FFTCircuit<Skew>::Apply(
            x0, x1, x2, x3, x4, x5, x6, x7,
            y0, y1, y2, y3, y4, y5, y6, y7);
    }

    ValueIO<ValueTag>::Store(x_buffer + offset, x0);
    ValueIO<ValueTag>::Store(x_buffer + plane_bytes + offset, x1);
    ValueIO<ValueTag>::Store(x_buffer + plane_bytes * 2 + offset, x2);
    ValueIO<ValueTag>::Store(x_buffer + plane_bytes * 3 + offset, x3);
    ValueIO<ValueTag>::Store(x_buffer + plane_bytes * 4 + offset, x4);
    ValueIO<ValueTag>::Store(x_buffer + plane_bytes * 5 + offset, x5);
    ValueIO<ValueTag>::Store(x_buffer + plane_bytes * 6 + offset, x6);
    ValueIO<ValueTag>::Store(x_buffer + plane_bytes * 7 + offset, x7);
    ValueIO<ValueTag>::Store(y_buffer + offset, y0);
    ValueIO<ValueTag>::Store(y_buffer + plane_bytes + offset, y1);
    ValueIO<ValueTag>::Store(y_buffer + plane_bytes * 2 + offset, y2);
    ValueIO<ValueTag>::Store(y_buffer + plane_bytes * 3 + offset, y3);
    ValueIO<ValueTag>::Store(y_buffer + plane_bytes * 4 + offset, y4);
    ValueIO<ValueTag>::Store(y_buffer + plane_bytes * 5 + offset, y5);
    ValueIO<ValueTag>::Store(y_buffer + plane_bytes * 6 + offset, y6);
    ValueIO<ValueTag>::Store(y_buffer + plane_bytes * 7 + offset, y7);
}


//------------------------------------------------------------------------------
// Coefficient-specialized whole-buffer kernels

template <unsigned Coefficient>
static void MultiplyWholeBuffer(
    void* destination_void,
    const void* source_void,
    uint64_t buffer_bytes)
{
    uint8_t* destination = reinterpret_cast<uint8_t*>(destination_void);
    const uint8_t* source = reinterpret_cast<const uint8_t*>(source_void);
    const uint64_t plane_bytes = buffer_bytes >> 3;
    const KernelMode mode = ResolveKernelMode();
    uint64_t offset = 0;

#ifdef LEO_TRY_AVX2
    if (mode == KernelMode::Avx2)
    {
        while (plane_bytes - offset >= ValueIO<Avx2Tag>::kBytes)
        {
            MultiplyChunk<Coefficient, Avx2Tag>(
                destination, source, plane_bytes, offset);
            offset += ValueIO<Avx2Tag>::kBytes;
        }
    }
#endif

#ifdef LEO_FF8XOR_HAS_SIMD128
    if (mode != KernelMode::Portable)
    {
        while (plane_bytes - offset >= ValueIO<Simd128Tag>::kBytes)
        {
            MultiplyChunk<Coefficient, Simd128Tag>(
                destination, source, plane_bytes, offset);
            offset += ValueIO<Simd128Tag>::kBytes;
        }
    }
#endif

    while (offset < plane_bytes)
    {
        MultiplyChunk<Coefficient, PortableTag>(
            destination, source, plane_bytes, offset);
        offset += ValueIO<PortableTag>::kBytes;
    }
}

template <unsigned Skew, bool Inverse>
static void ButterflyWholeBuffer(
    void* x_void,
    void* y_void,
    uint64_t buffer_bytes)
{
    uint8_t* x_buffer = reinterpret_cast<uint8_t*>(x_void);
    uint8_t* y_buffer = reinterpret_cast<uint8_t*>(y_void);
    const uint64_t plane_bytes = buffer_bytes >> 3;
    const KernelMode mode = ResolveKernelMode();
    uint64_t offset = 0;

#ifdef LEO_TRY_AVX2
    if (mode == KernelMode::Avx2)
    {
        while (plane_bytes - offset >= ValueIO<Avx2Tag>::kBytes)
        {
            ButterflyChunk<Skew, Avx2Tag, Inverse>(
                x_buffer, y_buffer, plane_bytes, offset);
            offset += ValueIO<Avx2Tag>::kBytes;
        }
    }
#endif

#ifdef LEO_FF8XOR_HAS_SIMD128
    if (mode != KernelMode::Portable)
    {
        while (plane_bytes - offset >= ValueIO<Simd128Tag>::kBytes)
        {
            ButterflyChunk<Skew, Simd128Tag, Inverse>(
                x_buffer, y_buffer, plane_bytes, offset);
            offset += ValueIO<Simd128Tag>::kBytes;
        }
    }
#endif

    while (offset < plane_bytes)
    {
        ButterflyChunk<Skew, PortableTag, Inverse>(
            x_buffer, y_buffer, plane_bytes, offset);
        offset += ValueIO<PortableTag>::kBytes;
    }
}

typedef void (*MultiplyFunction)(void*, const void*, uint64_t);
typedef void (*ButterflyFunction)(void*, void*, uint64_t);

#define LEO_FF8XOR_MULTIPLY_FUNCTION(Log) &MultiplyWholeBuffer<Log>,
static const MultiplyFunction MultiplyFunctions[kOrder] = {
    LEOPARD_FF8XOR_FOR_EACH_LOG(LEO_FF8XOR_MULTIPLY_FUNCTION)
};
#undef LEO_FF8XOR_MULTIPLY_FUNCTION

#define LEO_FF8XOR_FFT_FUNCTION(Skew) &ButterflyWholeBuffer<Skew, false>,
static const ButterflyFunction FFTFunctions[kOrder] = {
    LEOPARD_FF8XOR_FOR_EACH_LOG(LEO_FF8XOR_FFT_FUNCTION)
};
#undef LEO_FF8XOR_FFT_FUNCTION

#define LEO_FF8XOR_IFFT_FUNCTION(Skew) &ButterflyWholeBuffer<Skew, true>,
static const ButterflyFunction IFFTFunctions[kOrder] = {
    LEOPARD_FF8XOR_FOR_EACH_LOG(LEO_FF8XOR_IFFT_FUNCTION)
};
#undef LEO_FF8XOR_IFFT_FUNCTION

void MultiplyBuffer(
    uint64_t buffer_bytes,
    void* destination,
    const void* source,
    ffe_t log_multiplier)
{
    MultiplyFunctions[log_multiplier](destination, source, buffer_bytes);
}

void FFTButterflyBuffer(
    uint64_t buffer_bytes,
    void* x,
    void* y,
    ffe_t skew)
{
    FFTFunctions[skew](x, y, buffer_bytes);
}

void IFFTButterflyBuffer(
    uint64_t buffer_bytes,
    void* x,
    void* y,
    ffe_t skew)
{
    IFFTFunctions[skew](x, y, buffer_bytes);
}


//------------------------------------------------------------------------------
// Circuit metadata

const char* GetCircuitChecksum()
{
    return generated::kCircuitChecksum;
}

static CircuitStatistics MakeCircuitStatistics(
    unsigned minimum_gate_count,
    unsigned maximum_gate_count,
    double average_gate_count,
    unsigned minimum_depth,
    unsigned maximum_depth,
    double average_depth)
{
    CircuitStatistics statistics;
    statistics.MinimumGateCount = minimum_gate_count;
    statistics.MaximumGateCount = maximum_gate_count;
    statistics.AverageGateCount = average_gate_count;
    statistics.MinimumDepth = minimum_depth;
    statistics.MaximumDepth = maximum_depth;
    statistics.AverageDepth = average_depth;
    return statistics;
}

CircuitStatistics GetMultiplyCircuitStatistics()
{
    return MakeCircuitStatistics(
        generated::kMultiplyMinGateCount,
        generated::kMultiplyMaxGateCount,
        generated::kMultiplyAverageGateCount,
        generated::kMultiplyMinDepth,
        generated::kMultiplyMaxDepth,
        generated::kMultiplyAverageDepth);
}

CircuitStatistics GetFFTCircuitStatistics()
{
    return MakeCircuitStatistics(
        generated::kFFTMinGateCount,
        generated::kFFTMaxGateCount,
        generated::kFFTAverageGateCount,
        generated::kFFTMinDepth,
        generated::kFFTMaxDepth,
        generated::kFFTAverageDepth);
}

CircuitStatistics GetIFFTCircuitStatistics()
{
    return MakeCircuitStatistics(
        generated::kIFFTMinGateCount,
        generated::kIFFTMaxGateCount,
        generated::kIFFTAverageGateCount,
        generated::kIFFTMinDepth,
        generated::kIFFTMaxDepth,
        generated::kIFFTAverageDepth);
}


//------------------------------------------------------------------------------
// Scalar field metadata

static const ffe_t kCantorBasis[kBits] = {
    1, 214, 152, 146, 86, 200, 88, 230
};

static ffe_t LogLUT[kOrder];
static ffe_t ExpLUT[kOrder];

static LEO_FORCE_INLINE ffe_t AddMod(ffe_t a, ffe_t b)
{
    const unsigned sum = static_cast<unsigned>(a) + b;
    return static_cast<ffe_t>(sum + (sum >> kBits));
}

static LEO_FORCE_INLINE ffe_t SubMod(ffe_t a, ffe_t b)
{
    const unsigned difference = static_cast<unsigned>(a) - b;
    return static_cast<ffe_t>(difference + (difference >> kBits));
}

static ffe_t MultiplyLog(ffe_t value, ffe_t log_multiplier)
{
    if (value == 0)
        return 0;
    return ExpLUT[AddMod(LogLUT[value], log_multiplier)];
}

static void InitializeLogarithmTables()
{
    unsigned state = 1;
    for (unsigned exponent = 0; exponent < kModulus; ++exponent)
    {
        ExpLUT[state] = static_cast<ffe_t>(exponent);
        state <<= 1;
        if (state >= kOrder)
            state ^= kPolynomial;
    }
    ExpLUT[0] = kModulus;

    LogLUT[0] = 0;
    for (unsigned bit = 0; bit < kBits; ++bit)
    {
        const unsigned width = 1UL << bit;
        for (unsigned value = 0; value < width; ++value)
            LogLUT[value + width] = LogLUT[value] ^ kCantorBasis[bit];
    }

    for (unsigned value = 0; value < kOrder; ++value)
        LogLUT[value] = ExpLUT[LogLUT[value]];

    for (unsigned value = 0; value < kOrder; ++value)
        ExpLUT[LogLUT[value]] = static_cast<ffe_t>(value);

    ExpLUT[kModulus] = ExpLUT[0];
}


//------------------------------------------------------------------------------
// Fast Walsh-Hadamard transform for locator metadata

static LEO_FORCE_INLINE void FWHT_2(
    ffe_t& LEO_RESTRICT a,
    ffe_t& LEO_RESTRICT b)
{
    const ffe_t sum = AddMod(a, b);
    const ffe_t difference = SubMod(a, b);
    a = sum;
    b = difference;
}

static LEO_FORCE_INLINE void FWHT_4(ffe_t* data, unsigned stride)
{
    const unsigned stride2 = stride << 1;
    ffe_t x0 = data[0];
    ffe_t x1 = data[stride];
    ffe_t x2 = data[stride2];
    ffe_t x3 = data[stride2 + stride];

    FWHT_2(x0, x1);
    FWHT_2(x2, x3);
    FWHT_2(x0, x2);
    FWHT_2(x1, x3);

    data[0] = x0;
    data[stride] = x1;
    data[stride2] = x2;
    data[stride2 + stride] = x3;
}

static void FWHT(ffe_t* data, unsigned count, unsigned count_truncated)
{
    unsigned distance = 1;
    unsigned distance4 = 4;
    for (; distance4 <= count; distance = distance4, distance4 <<= 2)
    {
        for (unsigned range = 0; range < count_truncated; range += distance4)
        {
            for (unsigned index = range; index < range + distance; ++index)
                FWHT_4(data + index, distance);
        }
    }

    if (distance < count)
    {
        for (unsigned index = 0; index < distance; ++index)
            FWHT_2(data[index], data[index + distance]);
    }
}


//------------------------------------------------------------------------------
// FFT skew and locator tables

// Logical FFTSkew[i] is stored at FFTSkewPadded[i + 1].  Transform code can
// use a zero base for the old FFTSkew - 1 view without forming an invalid
// one-before-begin pointer.
static ffe_t FFTSkewPadded[kOrder];
static ffe_t LogWalsh[kOrder];

static LEO_FORCE_INLINE ffe_t& FFTSkewAt(unsigned index)
{
    return FFTSkewPadded[index + 1];
}

static void FFTInitialize()
{
    ffe_t temporary[kBits - 1];
    for (unsigned bit = 1; bit < kBits; ++bit)
        temporary[bit - 1] = static_cast<ffe_t>(1UL << bit);

    FFTSkewPadded[0] = kModulus;
    for (unsigned level = 0; level < kBits - 1; ++level)
    {
        const unsigned step = 1UL << (level + 1);
        FFTSkewAt((1UL << level) - 1) = 0;

        for (unsigned bit = level; bit < kBits - 1; ++bit)
        {
            const unsigned span = 1UL << (bit + 1);
            for (unsigned index = (1UL << level) - 1;
                index < span;
                index += step)
            {
                FFTSkewAt(index + span) = FFTSkewAt(index) ^ temporary[bit];
            }
        }

        temporary[level] = static_cast<ffe_t>(
            kModulus - LogLUT[MultiplyLog(
                temporary[level], LogLUT[temporary[level] ^ 1])]);

        for (unsigned bit = level + 1; bit < kBits - 1; ++bit)
        {
            const ffe_t sum = AddMod(
                LogLUT[temporary[bit] ^ 1], temporary[level]);
            temporary[bit] = MultiplyLog(temporary[bit], sum);
        }
    }

    for (unsigned index = 0; index < kModulus; ++index)
        FFTSkewAt(index) = LogLUT[FFTSkewAt(index)];

    for (unsigned index = 0; index < kOrder; ++index)
        LogWalsh[index] = LogLUT[index];
    LogWalsh[0] = 0;
    FWHT(LogWalsh, kOrder, kOrder);
}


//------------------------------------------------------------------------------
// Natural two-way transform schedules

static void IFFT_DIT(
    uint64_t buffer_bytes,
    unsigned count_truncated,
    void** work,
    unsigned count,
    unsigned skew_base)
{
    for (unsigned distance = 1; distance < count; distance <<= 1)
    {
        const unsigned span = distance << 1;
        for (unsigned range = 0;
            range < count_truncated;
            range += span)
        {
            const ffe_t skew = FFTSkewPadded[skew_base + range + distance];
            for (unsigned index = range; index < range + distance; ++index)
            {
                IFFTButterflyBuffer(
                    buffer_bytes,
                    work[index],
                    work[index + distance],
                    skew);
            }
        }
    }
}

static void FFT_DIT(
    uint64_t buffer_bytes,
    void** work,
    unsigned count_truncated,
    unsigned count,
    unsigned skew_base)
{
    for (unsigned distance = count >> 1; distance != 0; distance >>= 1)
    {
        const unsigned span = distance << 1;
        for (unsigned range = 0;
            range < count_truncated;
            range += span)
        {
            const ffe_t skew = FFTSkewPadded[skew_base + range + distance];
            for (unsigned index = range; index < range + distance; ++index)
            {
                FFTButterflyBuffer(
                    buffer_bytes,
                    work[index],
                    work[index + distance],
                    skew);
            }
        }
    }
}


//------------------------------------------------------------------------------
// Error-bitfield optimization

#ifdef LEO_ERROR_BITFIELD_OPT

class ErrorBitfield
{
    static const unsigned kWordCount = kOrder / 64;
    uint64_t Words[7][kWordCount];

public:
    ErrorBitfield()
    {
        memset(Words, 0, sizeof(Words));
    }

    LEO_FORCE_INLINE void Set(unsigned index)
    {
        Words[0][index / 64] |= static_cast<uint64_t>(1) << (index % 64);
    }

    void Prepare();

    LEO_FORCE_INLINE bool IsNeeded(unsigned mip_level, unsigned bit) const
    {
        if (mip_level >= 8)
            return true;
        return 0 != (Words[mip_level - 1][bit / 64]
            & (static_cast<uint64_t>(1) << (bit % 64)));
    }
};

static const uint64_t kHiMasks[5] = {
    0xAAAAAAAAAAAAAAAAULL,
    0xCCCCCCCCCCCCCCCCULL,
    0xF0F0F0F0F0F0F0F0ULL,
    0xFF00FF00FF00FF00ULL,
    0xFFFF0000FFFF0000ULL,
};

void ErrorBitfield::Prepare()
{
    for (unsigned word_index = 0; word_index < kWordCount; ++word_index)
    {
        uint64_t word = Words[0][word_index];
        const uint64_t high_to_low = word | ((word & kHiMasks[0]) >> 1);
        const uint64_t low_to_high = (word & (kHiMasks[0] >> 1)) << 1;
        Words[0][word_index] = word = high_to_low | low_to_high;

        for (unsigned level = 1, bits = 2; level < 5; ++level, bits <<= 1)
        {
            const uint64_t high = word | ((word & kHiMasks[level]) >> bits);
            const uint64_t low = (word & (kHiMasks[level] >> bits)) << bits;
            Words[level][word_index] = word = high | low;
        }
    }

    for (unsigned word_index = 0; word_index < kWordCount; ++word_index)
    {
        uint64_t word = Words[4][word_index];
        word |= word >> 32;
        word |= word << 32;
        Words[5][word_index] = word;
    }

    for (unsigned word_index = 0; word_index < kWordCount; word_index += 2)
    {
        const uint64_t word =
            Words[5][word_index] | Words[5][word_index + 1];
        Words[6][word_index] = word;
        Words[6][word_index + 1] = word;
    }
}

static void FFT_DIT_ErrorBits(
    uint64_t buffer_bytes,
    void** work,
    unsigned count_truncated,
    unsigned count,
    unsigned skew_base,
    const ErrorBitfield& error_bits)
{
    unsigned mip_level = LastNonzeroBit32(count);
    for (unsigned distance = count >> 1;
        distance != 0;
        distance >>= 1, --mip_level)
    {
        const unsigned span = distance << 1;
        for (unsigned range = 0;
            range < count_truncated;
            range += span)
        {
            if (!error_bits.IsNeeded(mip_level, range))
                continue;

            const ffe_t skew = FFTSkewPadded[skew_base + range + distance];
            for (unsigned index = range; index < range + distance; ++index)
            {
                FFTButterflyBuffer(
                    buffer_bytes,
                    work[index],
                    work[index + distance],
                    skew);
            }
        }
    }
}

#endif // LEO_ERROR_BITFIELD_OPT


//------------------------------------------------------------------------------
// Common error-locator log shift

static LEO_FORCE_INLINE unsigned CanonicalLog(ffe_t logarithm)
{
    return logarithm == kModulus ? 0 : logarithm;
}

static LEO_FORCE_INLINE ffe_t ShiftedLog(ffe_t logarithm, unsigned shift)
{
    unsigned result = CanonicalLog(logarithm) + shift;
    if (result >= kModulus)
        result -= kModulus;
    return static_cast<ffe_t>(result);
}

static LEO_FORCE_INLINE ffe_t InverseShiftedLog(
    ffe_t logarithm,
    unsigned shift)
{
    const unsigned shifted = ShiftedLog(logarithm, shift);
    return static_cast<ffe_t>(shifted == 0 ? 0 : kModulus - shifted);
}

static unsigned SelectLocatorShift(
    const ffe_t* error_locations,
    unsigned original_count,
    unsigned recovery_count,
    unsigned m,
    const void* const* original,
    const void* const* recovery)
{
    uint64_t best_gates = UINT64_MAX;
    uint64_t best_depth = UINT64_MAX;
    unsigned best_shift = 0;

    for (unsigned shift = 0; shift < kModulus; ++shift)
    {
        uint64_t gates = 0;
        uint64_t depth = 0;

        for (unsigned index = 0; index < recovery_count; ++index)
        {
            if (recovery[index])
            {
                const ffe_t log = ShiftedLog(error_locations[index], shift);
                gates += generated::kMultiplyGateCounts[log];
                depth += generated::kMultiplyDepths[log];
            }
        }

        for (unsigned index = 0; index < original_count; ++index)
        {
            if (original[index])
            {
                const ffe_t log = ShiftedLog(
                    error_locations[m + index], shift);
                gates += generated::kMultiplyGateCounts[log];
                depth += generated::kMultiplyDepths[log];
            }
            else
            {
                const ffe_t log = InverseShiftedLog(
                    error_locations[m + index], shift);
                gates += generated::kMultiplyGateCounts[log];
                depth += generated::kMultiplyDepths[log];
            }
        }

        if (gates < best_gates ||
            (gates == best_gates && depth < best_depth))
        {
            best_gates = gates;
            best_depth = depth;
            best_shift = shift;
        }
    }

    return best_shift;
}


//------------------------------------------------------------------------------
// Reed-Solomon encoding

void ReedSolomonEncode(
    uint64_t buffer_bytes,
    unsigned original_count,
    unsigned recovery_count,
    unsigned m,
    const void* const* data,
    void** work)
{
    unsigned chunk_start = 0;
    while (chunk_start < original_count)
    {
        const unsigned remaining = original_count - chunk_start;
        const unsigned chunk_count = remaining < m ? remaining : m;
        void** chunk_work = chunk_start == 0 ? work : work + m;

        for (unsigned index = 0; index < chunk_count; ++index)
        {
            memcpy(
                chunk_work[index],
                data[chunk_start + index],
                static_cast<size_t>(buffer_bytes));
        }
        for (unsigned index = chunk_count; index < m; ++index)
            memset(chunk_work[index], 0, static_cast<size_t>(buffer_bytes));

        // Padded index corresponding to logical FFTSkew + m - 1 + chunk_start.
        IFFT_DIT(
            buffer_bytes,
            chunk_count,
            chunk_work,
            m,
            m + chunk_start);

        if (chunk_start != 0)
            VectorXOR(buffer_bytes, m, work, chunk_work);

        chunk_start += m;
    }

    FFT_DIT(buffer_bytes, work, recovery_count, m, 0);
}


//------------------------------------------------------------------------------
// Reed-Solomon decoding

void ReedSolomonDecode(
    uint64_t buffer_bytes,
    unsigned original_count,
    unsigned recovery_count,
    unsigned m,
    unsigned n,
    const void* const* original,
    const void* const* recovery,
    void** work)
{
#ifdef LEO_ERROR_BITFIELD_OPT
    ErrorBitfield error_bits;
#endif

    ffe_t error_locations[kOrder] = {};
    for (unsigned index = 0; index < recovery_count; ++index)
    {
        if (!recovery[index])
            error_locations[index] = 1;
    }
    for (unsigned index = recovery_count; index < m; ++index)
        error_locations[index] = 1;

    for (unsigned index = 0; index < original_count; ++index)
    {
        if (!original[index])
        {
            error_locations[m + index] = 1;
#ifdef LEO_ERROR_BITFIELD_OPT
            error_bits.Set(m + index);
#endif
        }
    }

#ifdef LEO_ERROR_BITFIELD_OPT
    error_bits.Prepare();
#endif

    FWHT(error_locations, kOrder, m + original_count);
    for (unsigned index = 0; index < kOrder; ++index)
    {
        error_locations[index] = static_cast<ffe_t>(
            (static_cast<unsigned>(error_locations[index]) * LogWalsh[index])
            % kModulus);
    }
    FWHT(error_locations, kOrder, kOrder);

    const unsigned locator_shift = LocatorShiftOverride >= 0
        ? static_cast<unsigned>(LocatorShiftOverride)
        : SelectLocatorShift(
            error_locations,
            original_count,
            recovery_count,
            m,
            original,
            recovery);
    LastLocatorShift = locator_shift;

    for (unsigned index = 0; index < recovery_count; ++index)
    {
        if (recovery[index])
        {
            MultiplyBuffer(
                buffer_bytes,
                work[index],
                recovery[index],
                ShiftedLog(error_locations[index], locator_shift));
        }
        else
            memset(work[index], 0, static_cast<size_t>(buffer_bytes));
    }
    for (unsigned index = recovery_count; index < m; ++index)
        memset(work[index], 0, static_cast<size_t>(buffer_bytes));

    for (unsigned index = 0; index < original_count; ++index)
    {
        if (original[index])
        {
            MultiplyBuffer(
                buffer_bytes,
                work[m + index],
                original[index],
                ShiftedLog(error_locations[m + index], locator_shift));
        }
        else
            memset(work[m + index], 0, static_cast<size_t>(buffer_bytes));
    }
    for (unsigned index = m + original_count; index < n; ++index)
        memset(work[index], 0, static_cast<size_t>(buffer_bytes));

    IFFT_DIT(buffer_bytes, m + original_count, work, n, 0);

    for (unsigned index = 1; index < n; ++index)
    {
        const unsigned width = ((index ^ (index - 1)) + 1) >> 1;
        VectorXOR(
            buffer_bytes,
            width,
            work + index - width,
            work + index);
    }

    const unsigned output_count = m + original_count;
#ifdef LEO_ERROR_BITFIELD_OPT
    FFT_DIT_ErrorBits(
        buffer_bytes,
        work,
        output_count,
        n,
        0,
        error_bits);
#else
    FFT_DIT(buffer_bytes, work, output_count, n, 0);
#endif

    for (unsigned index = 0; index < original_count; ++index)
    {
        if (!original[index])
        {
            MultiplyBuffer(
                buffer_bytes,
                work[index],
                work[index + m],
                InverseShiftedLog(
                    error_locations[index + m], locator_shift));
        }
    }
}


//------------------------------------------------------------------------------
// Initialization

static bool Initialized = false;

bool Initialize()
{
    if (Initialized)
        return true;

    InitializeLogarithmTables();
    FFTInitialize();
    Initialized = true;
    return true;
}

bool IsInitialized()
{
    return Initialized;
}


}} // namespace leopard::ff8xor

#endif // LEO_HAS_FF8
