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
#include "LeopardFF8XorAVX2.h"
#include "LeopardFF8XorAVX512.h"
#include "LeopardFF8XorAVX512Four.h"
#include "LeopardFF8XorDerivative.h"

#ifdef LEO_HAS_FF8

#include <atomic>
#include <limits.h>
#include <string.h>

#if (defined(LEO_FF8XOR_HAS_AVX2_KERNELS) || \
     defined(LEO_FF8XOR_HAS_AVX512_KERNELS)) && \
    defined(_MSC_VER) && \
    (defined(_M_X64) || defined(_M_AMD64) || defined(_M_IX86))
    #include <intrin.h>
#elif (defined(LEO_FF8XOR_HAS_AVX2_KERNELS) || \
       defined(LEO_FF8XOR_HAS_AVX512_KERNELS)) && \
    (defined(__i386__) || defined(__x86_64__))
    #include <cpuid.h>
#endif

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

}} // namespace leopard::ff8xor

#include "generated/LeopardFF8XorCircuits.inl"
#include "generated/LeopardFF8XorLocatorRotations.inl"

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

template <typename ValueTag>
struct BaselineDerivativeOps
{
    typedef typename ValueIO<ValueTag>::Value Value;
    static const unsigned kBytes = ValueIO<ValueTag>::kBytes;

    static LEO_FORCE_INLINE Value Load(const uint8_t* source)
    {
        return ValueIO<ValueTag>::Load(source);
    }

    static LEO_FORCE_INLINE void Store(uint8_t* destination, Value value)
    {
        ValueIO<ValueTag>::Store(destination, value);
    }

    static LEO_FORCE_INLINE Value Xor(Value a, Value b)
    {
        return XorValue(a, b);
    }

    static LEO_FORCE_INLINE Value Xor3(Value a, Value b, Value c)
    {
        return XorValue(XorValue(a, b), c);
    }
};

//------------------------------------------------------------------------------
// Kernel selection

// Test inspection and overrides must not introduce shared mutable state into
// otherwise independent codec calls.
static std::atomic<KernelMode> SelectedKernelMode(KernelMode::Auto);
static std::atomic<FourBufferMode> SelectedFourBufferMode(
    FourBufferMode::Disabled);
static thread_local unsigned LastLocatorShift = 0;
static thread_local int LocatorShiftOverride = -1;
static thread_local MaterializationStatistics LastMaterializationStatistics = {};
static thread_local FourBufferStatistics LastFourBufferStatistics = {};

bool IsAVX512Supported(
    uint32_t maximum_basic_leaf,
    uint32_t leaf1_ecx,
    uint32_t leaf7_ebx,
    uint64_t xcr0)
{
    const uint32_t required_leaf1 = UINT32_C(1) << 27 | // OSXSAVE
        UINT32_C(1) << 28;                              // AVX
    const uint32_t required_leaf7 = UINT32_C(1) << 5 |  // AVX2 tails
        UINT32_C(1) << 16 |                             // AVX-512F
        UINT32_C(1) << 31;                              // AVX-512VL
    return maximum_basic_leaf >= 7 &&
        (leaf1_ecx & required_leaf1) == required_leaf1 &&
        (leaf7_ebx & required_leaf7) == required_leaf7 &&
        (xcr0 & UINT64_C(0xe6)) == UINT64_C(0xe6);
}

#if defined(LEO_FF8XOR_HAS_AVX2_KERNELS) || \
    defined(LEO_FF8XOR_HAS_AVX512_KERNELS)

struct X86VectorCapabilities
{
    bool Avx2;
    bool Avx512Foundation;
    bool Avx512VectorLength;
};

// This detector is deliberately compiled in the baseline translation unit.
// No AVX-targeted object is entered until CPUID and XCR0 prove that both the
// processor and operating system preserve the vector state that object needs.
#if defined(_MSC_VER)
    #define LEO_FF8XOR_NOINLINE __declspec(noinline)
#elif defined(__GNUC__) || defined(__clang__)
    #define LEO_FF8XOR_NOINLINE __attribute__((noinline))
#else
    #define LEO_FF8XOR_NOINLINE
#endif

static LEO_FF8XOR_NOINLINE X86VectorCapabilities
DetectX86VectorCapabilities()
{
    X86VectorCapabilities capabilities = { false, false, false };

#if defined(_MSC_VER) && \
    (defined(_M_X64) || defined(_M_AMD64) || defined(_M_IX86))
    int registers[4] = { 0, 0, 0, 0 };
    __cpuidex(registers, 0, 0);
    const uint32_t maximum_basic_leaf =
        static_cast<uint32_t>(registers[0]);
    if (maximum_basic_leaf < 7)
        return capabilities;

    __cpuidex(registers, 1, 0);
    const unsigned leaf1_ecx = static_cast<unsigned>(registers[2]);
    if ((leaf1_ecx & (1U << 27)) == 0 || // OSXSAVE
        (leaf1_ecx & (1U << 28)) == 0)   // AVX
        return capabilities;

    __cpuidex(registers, 7, 0);
    const unsigned leaf7_ebx = static_cast<unsigned>(registers[1]);
    const uint64_t xcr0 = static_cast<uint64_t>(_xgetbv(0));
    capabilities.Avx2 = IsAVX2Supported(
        maximum_basic_leaf,
        leaf1_ecx,
        leaf7_ebx,
        xcr0);
    if (IsAVX512Supported(
            maximum_basic_leaf,
            leaf1_ecx,
            leaf7_ebx,
            xcr0))
    {
        capabilities.Avx512Foundation =
            (leaf7_ebx & (1U << 16)) != 0;
        capabilities.Avx512VectorLength =
            (leaf7_ebx & (1U << 31)) != 0;
    }
#elif defined(__i386__) || defined(__x86_64__)
    const uint32_t maximum_basic_leaf = __get_cpuid_max(0, NULL);
    if (maximum_basic_leaf < 7)
        return capabilities;

    unsigned eax = 0;
    unsigned ebx = 0;
    unsigned ecx = 0;
    unsigned edx = 0;
    __cpuid_count(1, 0, eax, ebx, ecx, edx);
    const uint32_t leaf1_ecx = ecx;
    if ((leaf1_ecx & (1U << 27)) == 0 || // OSXSAVE
        (leaf1_ecx & (1U << 28)) == 0)   // AVX
        return capabilities;

    unsigned xcr0_low = 0;
    unsigned xcr0_high = 0;
    // Raw encoding avoids requiring a compiler-wide -mxsave option.  This is
    // reached only after CPUID reports both AVX and OSXSAVE.
    __asm__ __volatile__(
        ".byte 0x0f, 0x01, 0xd0"
        : "=a"(xcr0_low), "=d"(xcr0_high) : "c"(0));
    const uint64_t xcr0 = static_cast<uint64_t>(xcr0_low) |
        (static_cast<uint64_t>(xcr0_high) << 32);
    __cpuid_count(7, 0, eax, ebx, ecx, edx);
    capabilities.Avx2 = IsAVX2Supported(
        maximum_basic_leaf, leaf1_ecx, ebx, xcr0);
    if (IsAVX512Supported(
            maximum_basic_leaf, leaf1_ecx, ebx, xcr0))
    {
        capabilities.Avx512Foundation = (ebx & (1U << 16)) != 0;
        capabilities.Avx512VectorLength = (ebx & (1U << 31)) != 0;
    }
#endif

    return capabilities;
}

#undef LEO_FF8XOR_NOINLINE

static const X86VectorCapabilities& GetX86VectorCapabilities()
{
    static const X86VectorCapabilities capabilities =
        DetectX86VectorCapabilities();
    return capabilities;
}

#endif

static bool AreAVX2KernelsAvailable()
{
#ifdef LEO_FF8XOR_HAS_AVX2_KERNELS
    // Cache the immutable result so whole-buffer coefficient dispatch never
    // pays an extra cross-TU feature-query call.
    static const bool available =
        GetX86VectorCapabilities().Avx2 && avx2::KernelsBuilt();
    return available;
#else
    return false;
#endif
}

#ifdef LEO_FF8XOR_HAS_AVX512_KERNELS

static bool DetectAVX512KernelAvailability()
{
    const X86VectorCapabilities& capabilities =
        GetX86VectorCapabilities();

    // The isolated TU is compiled with AVX2, F, and VL enabled.  Requiring all
    // three features for either width is conservative: It prevents a compiler
    // from using a permitted AVX2/VL instruction in an otherwise ZMM-only
    // function on a hypothetical AVX-512F CPU lacking either feature.
    return AreAVX2KernelsAvailable() && capabilities.Avx512Foundation &&
        capabilities.Avx512VectorLength && avx512::KernelsBuilt();
}

#endif // LEO_FF8XOR_HAS_AVX512_KERNELS

static bool AreAVX512KernelsAvailable()
{
#ifdef LEO_FF8XOR_HAS_AVX512_KERNELS
    // Dispatch is per whole-buffer operation.  Cache the immutable result so
    // it does not add a cross-TU availability call to every coefficient.
    static const bool available = DetectAVX512KernelAvailability();
    return available;
#else
    return false;
#endif
}

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
        return AreAVX2KernelsAvailable();

    case KernelMode::Avx512VL:
    case KernelMode::Avx512Zmm:
        return AreAVX512KernelsAvailable();
    }
    return false;
}

static KernelMode ResolveBaselineKernelMode()
{
    if (AreAVX2KernelsAvailable())
        return KernelMode::Avx2;
#ifdef LEO_FF8XOR_HAS_SIMD128
    return KernelMode::Simd128;
#else
    return KernelMode::Portable;
#endif
}

static KernelMode ResolveKernelMode()
{
    const KernelMode selected =
        SelectedKernelMode.load(std::memory_order_relaxed);

    if (selected == KernelMode::Portable)
        return KernelMode::Portable;

    if (selected == KernelMode::Avx512VL ||
        selected == KernelMode::Avx512Zmm)
    {
        if (AreAVX512KernelsAvailable())
            return selected;
        return ResolveBaselineKernelMode();
    }

    if (selected == KernelMode::Avx2)
    {
        if (AreAVX2KernelsAvailable())
            return KernelMode::Avx2;
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

    // Keep Auto on the established backend until comparative complete-codec
    // measurements select an AVX-512 mode without regressing other workloads.
    return ResolveBaselineKernelMode();
}

void SetKernelMode(KernelMode mode)
{
    SelectedKernelMode.store(mode, std::memory_order_relaxed);
}

void SetFourBufferMode(FourBufferMode mode)
{
    SelectedFourBufferMode.store(mode, std::memory_order_relaxed);
}

FourBufferMode GetFourBufferMode()
{
    return SelectedFourBufferMode.load(std::memory_order_relaxed);
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
    case KernelMode::Avx512Zmm: return "AVX-512 ZMM XOR circuits";
    case KernelMode::Avx512VL: return "AVX-512VL YMM XOR circuits";
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

MaterializationStatistics GetLastMaterializationStatistics()
{
    return LastMaterializationStatistics;
}

void ResetMaterializationStatistics()
{
    const MaterializationStatistics empty = {};
    LastMaterializationStatistics = empty;
}

FourBufferStatistics GetLastFourBufferStatistics()
{
    return LastFourBufferStatistics;
}

void ResetFourBufferStatistics()
{
    const FourBufferStatistics empty = {};
    LastFourBufferStatistics = empty;
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

template <typename ValueTag>
static LEO_FORCE_INLINE void XorContiguousChunk(
    uint8_t* destination,
    const uint8_t* source,
    uint64_t offset)
{
    typedef typename ValueIO<ValueTag>::Value Value;
    const Value result = XorValue(
        ValueIO<ValueTag>::Load(destination + offset),
        ValueIO<ValueTag>::Load(source + offset));
    ValueIO<ValueTag>::Store(destination + offset, result);
}

template <typename ValueTag>
static LEO_FORCE_INLINE void XorContiguous2Chunk(
    uint8_t* destination0,
    const uint8_t* source0,
    uint8_t* destination1,
    const uint8_t* source1,
    uint64_t offset)
{
    typedef typename ValueIO<ValueTag>::Value Value;
    const Value result0 = XorValue(
        ValueIO<ValueTag>::Load(destination0 + offset),
        ValueIO<ValueTag>::Load(source0 + offset));
    const Value result1 = XorValue(
        ValueIO<ValueTag>::Load(destination1 + offset),
        ValueIO<ValueTag>::Load(source1 + offset));
    ValueIO<ValueTag>::Store(destination0 + offset, result0);
    ValueIO<ValueTag>::Store(destination1 + offset, result1);
}

template <typename ValueTag>
static LEO_FORCE_INLINE void XorContiguous4Chunk(
    uint8_t* destination0,
    const uint8_t* source0,
    uint8_t* destination1,
    const uint8_t* source1,
    uint8_t* destination2,
    const uint8_t* source2,
    uint8_t* destination3,
    const uint8_t* source3,
    uint64_t offset)
{
    typedef typename ValueIO<ValueTag>::Value Value;
    const Value result0 = XorValue(
        ValueIO<ValueTag>::Load(destination0 + offset),
        ValueIO<ValueTag>::Load(source0 + offset));
    const Value result1 = XorValue(
        ValueIO<ValueTag>::Load(destination1 + offset),
        ValueIO<ValueTag>::Load(source1 + offset));
    const Value result2 = XorValue(
        ValueIO<ValueTag>::Load(destination2 + offset),
        ValueIO<ValueTag>::Load(source2 + offset));
    const Value result3 = XorValue(
        ValueIO<ValueTag>::Load(destination3 + offset),
        ValueIO<ValueTag>::Load(source3 + offset));
    ValueIO<ValueTag>::Store(destination0 + offset, result0);
    ValueIO<ValueTag>::Store(destination1 + offset, result1);
    ValueIO<ValueTag>::Store(destination2 + offset, result2);
    ValueIO<ValueTag>::Store(destination3 + offset, result3);
}

// Skew 255 is not a multiplier.  It means that the multiply-add term is
// omitted, leaving x unchanged and y ^= x.  Process that operation as one
// contiguous buffer so the payload loop never stores x and does not pay the
// eight-plane address-arithmetic overhead of a general butterfly.
static void XorContiguousWholeBuffer(
    KernelMode mode,
    void* destination_void,
    const void* source_void,
    uint64_t buffer_bytes)
{
    uint8_t* destination = reinterpret_cast<uint8_t*>(destination_void);
    const uint8_t* source = reinterpret_cast<const uint8_t*>(source_void);
    uint64_t offset = 0;

#ifdef LEO_FF8XOR_HAS_AVX512_KERNELS
    if (mode == KernelMode::Avx512Zmm)
        offset = avx512::Xor512(destination, source, buffer_bytes);
    else if (mode == KernelMode::Avx512VL)
        offset = avx512::Xor256(destination, source, buffer_bytes);
#endif

#ifdef LEO_FF8XOR_HAS_AVX2_KERNELS
    if (mode == KernelMode::Avx2 ||
        mode == KernelMode::Avx512VL ||
        mode == KernelMode::Avx512Zmm)
        offset = avx2::Xor(
            destination, source, buffer_bytes, offset);
#endif

#ifdef LEO_FF8XOR_HAS_SIMD128
    if (mode != KernelMode::Portable)
    {
        while (buffer_bytes - offset >= ValueIO<Simd128Tag>::kBytes)
        {
            XorContiguousChunk<Simd128Tag>(destination, source, offset);
            offset += ValueIO<Simd128Tag>::kBytes;
        }
    }
#endif

    while (offset < buffer_bytes)
    {
        XorContiguousChunk<PortableTag>(destination, source, offset);
        offset += ValueIO<PortableTag>::kBytes;
    }
}

static void XorContiguousWholeBuffer(
    void* destination,
    const void* source,
    uint64_t buffer_bytes)
{
    XorContiguousWholeBuffer(
        ResolveKernelMode(), destination, source, buffer_bytes);
}

static void XorContiguous2WholeBuffers(
    KernelMode mode,
    void* destination0_void,
    const void* source0_void,
    void* destination1_void,
    const void* source1_void,
    uint64_t buffer_bytes)
{
    uint8_t* destination0 = reinterpret_cast<uint8_t*>(destination0_void);
    const uint8_t* source0 = reinterpret_cast<const uint8_t*>(source0_void);
    uint8_t* destination1 = reinterpret_cast<uint8_t*>(destination1_void);
    const uint8_t* source1 = reinterpret_cast<const uint8_t*>(source1_void);
    uint64_t offset = 0;

#ifdef LEO_FF8XOR_HAS_AVX512_KERNELS
    if (mode == KernelMode::Avx512Zmm)
        offset = avx512::Xor2_512(
            destination0, source0, destination1, source1, buffer_bytes);
    else if (mode == KernelMode::Avx512VL)
        offset = avx512::Xor2_256(
            destination0, source0, destination1, source1, buffer_bytes);
#endif

#ifdef LEO_FF8XOR_HAS_AVX2_KERNELS
    if (mode == KernelMode::Avx2 ||
        mode == KernelMode::Avx512VL ||
        mode == KernelMode::Avx512Zmm)
    {
        offset = avx2::Xor2(
            destination0, source0, destination1, source1,
            buffer_bytes, offset);
    }
#endif

#ifdef LEO_FF8XOR_HAS_SIMD128
    if (mode != KernelMode::Portable)
    {
        while (buffer_bytes - offset >= ValueIO<Simd128Tag>::kBytes)
        {
            XorContiguous2Chunk<Simd128Tag>(
                destination0, source0, destination1, source1, offset);
            offset += ValueIO<Simd128Tag>::kBytes;
        }
    }
#endif

    while (offset < buffer_bytes)
    {
        XorContiguous2Chunk<PortableTag>(
            destination0, source0, destination1, source1, offset);
        offset += ValueIO<PortableTag>::kBytes;
    }
}

static void XorContiguous4WholeBuffers(
    KernelMode mode,
    void* destination0_void,
    const void* source0_void,
    void* destination1_void,
    const void* source1_void,
    void* destination2_void,
    const void* source2_void,
    void* destination3_void,
    const void* source3_void,
    uint64_t buffer_bytes)
{
    uint8_t* destination0 = reinterpret_cast<uint8_t*>(destination0_void);
    const uint8_t* source0 = reinterpret_cast<const uint8_t*>(source0_void);
    uint8_t* destination1 = reinterpret_cast<uint8_t*>(destination1_void);
    const uint8_t* source1 = reinterpret_cast<const uint8_t*>(source1_void);
    uint8_t* destination2 = reinterpret_cast<uint8_t*>(destination2_void);
    const uint8_t* source2 = reinterpret_cast<const uint8_t*>(source2_void);
    uint8_t* destination3 = reinterpret_cast<uint8_t*>(destination3_void);
    const uint8_t* source3 = reinterpret_cast<const uint8_t*>(source3_void);
    uint64_t offset = 0;

#ifdef LEO_FF8XOR_HAS_AVX512_KERNELS
    if (mode == KernelMode::Avx512Zmm)
    {
        offset = avx512::Xor4_512(
            destination0, source0, destination1, source1,
            destination2, source2, destination3, source3, buffer_bytes);
    }
    else if (mode == KernelMode::Avx512VL)
    {
        offset = avx512::Xor4_256(
            destination0, source0, destination1, source1,
            destination2, source2, destination3, source3, buffer_bytes);
    }
#endif

#ifdef LEO_FF8XOR_HAS_AVX2_KERNELS
    if (mode == KernelMode::Avx2 ||
        mode == KernelMode::Avx512VL ||
        mode == KernelMode::Avx512Zmm)
    {
        offset = avx2::Xor4(
            destination0, source0, destination1, source1,
            destination2, source2, destination3, source3,
            buffer_bytes, offset);
    }
#endif

#ifdef LEO_FF8XOR_HAS_SIMD128
    if (mode != KernelMode::Portable)
    {
        while (buffer_bytes - offset >= ValueIO<Simd128Tag>::kBytes)
        {
            XorContiguous4Chunk<Simd128Tag>(
                destination0, source0, destination1, source1,
                destination2, source2, destination3, source3, offset);
            offset += ValueIO<Simd128Tag>::kBytes;
        }
    }
#endif

    while (offset < buffer_bytes)
    {
        XorContiguous4Chunk<PortableTag>(
            destination0, source0, destination1, source1,
            destination2, source2, destination3, source3, offset);
        offset += ValueIO<PortableTag>::kBytes;
    }
}

void XorBuffer(
    uint64_t buffer_bytes,
    void* destination,
    const void* source)
{
    XorContiguousWholeBuffer(destination, source, buffer_bytes);
}

void XorBuffers(
    uint64_t buffer_bytes,
    unsigned count,
    void** destination,
    void** source)
{
    if (count == 0)
        return;

    const KernelMode mode = ResolveKernelMode();
    unsigned index = 0;
    for (; count - index >= 4; index += 4)
    {
        if (buffer_bytes <= 1024)
        {
            XorContiguous4WholeBuffers(
                mode,
                destination[index], source[index],
                destination[index + 1], source[index + 1],
                destination[index + 2], source[index + 2],
                destination[index + 3], source[index + 3],
                buffer_bytes);
        }
        else
        {
            XorContiguousWholeBuffer(
                mode, destination[index], source[index], buffer_bytes);
            XorContiguousWholeBuffer(
                mode,
                destination[index + 1], source[index + 1], buffer_bytes);
            XorContiguousWholeBuffer(
                mode,
                destination[index + 2], source[index + 2], buffer_bytes);
            XorContiguousWholeBuffer(
                mode,
                destination[index + 3], source[index + 3], buffer_bytes);
        }
    }
    if (count - index >= 2)
    {
        // The two-stream loop removes dispatch overhead for small buffers, but
        // on this implementation host it slightly reduced sustained bandwidth
        // once each stream exceeded 1 KiB.  Preserve the single mode lookup
        // while using the independently tuned loop shape in each region.
        if (buffer_bytes <= 1024)
        {
            XorContiguous2WholeBuffers(
                mode,
                destination[index], source[index],
                destination[index + 1], source[index + 1],
                buffer_bytes);
        }
        else
        {
            XorContiguousWholeBuffer(
                mode, destination[index], source[index], buffer_bytes);
            XorContiguousWholeBuffer(
                mode,
                destination[index + 1], source[index + 1], buffer_bytes);
        }
        index += 2;
    }
    if (index < count)
    {
        XorContiguousWholeBuffer(
            mode, destination[index], source[index], buffer_bytes);
    }
}

static void ApplyFormalDerivativeBoundaryRow(
    KernelMode mode,
    unsigned extra_count,
    void* left,
    void* right,
    const void* extra0,
    const void* extra1,
    const void* extra2,
    const void* extra3,
    const void* extra4,
    const void* extra5,
    const void* extra6,
    uint64_t buffer_bytes)
{
    uint64_t offset = 0;

#ifdef LEO_FF8XOR_HAS_AVX512_KERNELS
    if (mode == KernelMode::Avx512Zmm)
    {
        offset = avx512::FormalDerivativeBoundaryRow512(
            extra_count, left, right,
            extra0, extra1, extra2, extra3, extra4, extra5, extra6,
            buffer_bytes);
    }
    else if (mode == KernelMode::Avx512VL)
    {
        offset = avx512::FormalDerivativeBoundaryRow256(
            extra_count, left, right,
            extra0, extra1, extra2, extra3, extra4, extra5, extra6,
            buffer_bytes);
    }
#endif

#ifdef LEO_FF8XOR_HAS_AVX2_KERNELS
    if (mode == KernelMode::Avx2 ||
        mode == KernelMode::Avx512VL ||
        mode == KernelMode::Avx512Zmm)
    {
        offset = avx2::FormalDerivativeBoundaryRow(
            extra_count, left, right,
            extra0, extra1, extra2, extra3, extra4, extra5, extra6,
            buffer_bytes, offset);
    }
#endif

#ifdef LEO_FF8XOR_HAS_SIMD128
    if (mode != KernelMode::Portable)
    {
        offset = derivative_detail::ApplyRow<
            BaselineDerivativeOps<Simd128Tag> >(
                extra_count, left, right,
                extra0, extra1, extra2, extra3, extra4, extra5, extra6,
                buffer_bytes, offset);
    }
#endif

    derivative_detail::ApplyRow<BaselineDerivativeOps<PortableTag> >(
        extra_count, left, right,
        extra0, extra1, extra2, extra3, extra4, extra5, extra6,
        buffer_bytes, offset);
}

// A is the already-computed derivative of the left half and R is the original
// right half.  The old schedule next computed D(R), XORed R into A, then ran
// the top FFT butterfly.  Its skew is sentinel 255 for every supported power
// of two.  The exact combined map is:
//
//   left[q]  = A[q] ^ R[q]
//   right[q] = A[q] ^ XOR(R[q | (1 << b)]) for every zero bit b of q
//
// Higher right-half indices are the only extra sources, so ascending q is
// safely in place.  Arity dispatch and pointer selection happen once per
// whole-buffer row, outside every payload loop.
static void ApplyFormalDerivativeTopFFTBoundary(
    uint64_t buffer_bytes,
    unsigned half_count,
    void** work)
{
    const KernelMode mode = ResolveKernelMode();
    for (unsigned q = 0; q < half_count; ++q)
    {
        const void* extra0 = NULL;
        const void* extra1 = NULL;
        const void* extra2 = NULL;
        const void* extra3 = NULL;
        const void* extra4 = NULL;
        const void* extra5 = NULL;
        const void* extra6 = NULL;
        unsigned extra_count = 0;

        for (unsigned bit = 1; bit < half_count; bit <<= 1)
        {
            if ((q & bit) != 0)
                continue;
            const void* source = work[half_count + (q | bit)];
            switch (extra_count++)
            {
            case 0: extra0 = source; break;
            case 1: extra1 = source; break;
            case 2: extra2 = source; break;
            case 3: extra3 = source; break;
            case 4: extra4 = source; break;
            case 5: extra5 = source; break;
            case 6: extra6 = source; break;
            default: break;
            }
        }

        ApplyFormalDerivativeBoundaryRow(
            mode,
            extra_count,
            work[q],
            work[half_count + q],
            extra0,
            extra1,
            extra2,
            extra3,
            extra4,
            extra5,
            extra6,
            buffer_bytes);
    }
}

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

#ifdef LEO_FF8XOR_HAS_AVX512_KERNELS
    if (mode == KernelMode::Avx512Zmm)
    {
        offset = avx512::Multiply512<Coefficient>(
            destination, source, buffer_bytes);
    }
    else if (mode == KernelMode::Avx512VL)
    {
        offset = avx512::Multiply256<Coefficient>(
            destination, source, buffer_bytes);
    }
#endif

#ifdef LEO_FF8XOR_HAS_AVX2_KERNELS
    if (mode == KernelMode::Avx2 ||
        mode == KernelMode::Avx512VL ||
        mode == KernelMode::Avx512Zmm)
        offset = avx2::Multiply<Coefficient>(
            destination, source, buffer_bytes, offset);
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
    // This constant branch is eliminated separately for every specialization.
    // Both FFT and IFFT have the same skew-sentinel operation.
    if (Skew == kModulus)
    {
        XorContiguousWholeBuffer(y_void, x_void, buffer_bytes);
        return;
    }

    uint8_t* x_buffer = reinterpret_cast<uint8_t*>(x_void);
    uint8_t* y_buffer = reinterpret_cast<uint8_t*>(y_void);
    const uint64_t plane_bytes = buffer_bytes >> 3;
    const KernelMode mode = ResolveKernelMode();
    uint64_t offset = 0;

#ifdef LEO_FF8XOR_HAS_AVX512_KERNELS
    if (mode == KernelMode::Avx512Zmm)
    {
        offset = avx512::Butterfly512<Skew, Inverse>(
            x_buffer, y_buffer, buffer_bytes);
    }
    else if (mode == KernelMode::Avx512VL)
    {
        offset = avx512::Butterfly256<Skew, Inverse>(
            x_buffer, y_buffer, buffer_bytes);
    }
#endif

#ifdef LEO_FF8XOR_HAS_AVX2_KERNELS
    if (mode == KernelMode::Avx2 ||
        mode == KernelMode::Avx512VL ||
        mode == KernelMode::Avx512Zmm)
        offset = avx2::Butterfly<Skew, Inverse>(
            x_buffer, y_buffer, buffer_bytes, offset);
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
    // Both logarithms denote the nonzero field element one for multiplication.
    // Keep this distinct from butterfly skew 255, which omits a multiply-add.
    // Exact in-place identity must not access payload memory at all.  memmove
    // gives the out-of-place identity operation defined copy semantics even
    // when the source and destination ranges overlap.
    if (log_multiplier == 0 || log_multiplier == kModulus)
    {
        if (destination != source)
            memmove(destination, source, static_cast<size_t>(buffer_bytes));
        return;
    }

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

const char* GetCircuitCostProfileId()
{
    return generated::kCircuitCostProfileId;
}

const char* GetCircuitCostProfileChecksum()
{
    return generated::kCircuitCostProfileChecksum;
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

static CircuitCost MakeCircuitCost(unsigned gate_count, unsigned depth)
{
    CircuitCost cost;
    cost.GateCount = gate_count;
    cost.Depth = depth;
    return cost;
}

CircuitCost GetMultiplyCircuitCost(ffe_t log_multiplier)
{
    return MakeCircuitCost(
        generated::kMultiplyGateCounts[log_multiplier],
        generated::kMultiplyDepths[log_multiplier]);
}

CircuitCost GetFFTCircuitCost(ffe_t skew)
{
    return MakeCircuitCost(
        generated::kFFTGateCounts[skew],
        generated::kFFTDepths[skew]);
}

CircuitCost GetIFFTCircuitCost(ffe_t skew)
{
    return MakeCircuitCost(
        generated::kIFFTGateCounts[skew],
        generated::kIFFTDepths[skew]);
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

struct BufferState
{
    // Zero is the distinguished logical-zero identity.  Nonzero identities
    // prove exact equality only when produced by a copy; arbitrary inputs get
    // distinct identities even if their payload bytes happen to match.
    uint32_t Identity;
    bool Materialized;
};

struct MaterializationContext
{
    uint64_t BufferBytes;
    uint32_t NextIdentity;
    MaterializationStatistics Statistics;

    explicit MaterializationContext(uint64_t buffer_bytes)
        : BufferBytes(buffer_bytes)
        , NextIdentity(1)
        , Statistics()
    {
    }

    void SetFresh(BufferState& state)
    {
        state.Identity = NextIdentity++;
        state.Materialized = true;
    }

    void SetZero(BufferState& state, bool materialized)
    {
        state.Identity = 0;
        state.Materialized = materialized;
    }

    void DeferZero(BufferState& state)
    {
        SetZero(state, false);
        ++Statistics.DeferredZeroFills;
        Statistics.EstimatedPayloadBytesElided +=
            static_cast<int64_t>(BufferBytes);
    }

    void MaterializeZero(void* buffer, BufferState& state)
    {
        if (state.Identity != 0 || state.Materialized)
            return;
        memset(buffer, 0, static_cast<size_t>(BufferBytes));
        state.Materialized = true;
        ++Statistics.AddedZeroFills;
        Statistics.EstimatedPayloadBytesElided -=
            static_cast<int64_t>(BufferBytes);
    }

    uint64_t ButterflyPayloadBytes(ffe_t skew) const
    {
        // Sentinel 255 is the existing contiguous y ^= x path: two reads and
        // one store.  A generated non-sentinel butterfly reads and stores both
        // buffers.
        return BufferBytes * (skew == kModulus ? 3 : 4);
    }

    void RecordSkippedButterfly(ffe_t skew, bool identity)
    {
        ++Statistics.ButterfliesSkipped;
        if (identity)
            ++Statistics.IdentityOperationsElided;
        Statistics.EstimatedPayloadBytesElided +=
            static_cast<int64_t>(ButterflyPayloadBytes(skew));
    }

    void RecordReducedButterfly(ffe_t skew, uint64_t replacement_bytes)
    {
        ++Statistics.ButterfliesReduced;
        Statistics.EstimatedPayloadBytesElided += static_cast<int64_t>(
            ButterflyPayloadBytes(skew) - replacement_bytes);
    }

    void RecordSkippedXor(bool identity)
    {
        ++Statistics.XorsSkipped;
        if (identity)
            ++Statistics.IdentityOperationsElided;
        Statistics.EstimatedPayloadBytesElided +=
            static_cast<int64_t>(BufferBytes * 3);
    }
};

static LEO_FORCE_INLINE bool IsZero(const BufferState& state)
{
    return state.Identity == 0;
}

static LEO_FORCE_INLINE void CopyBufferState(
    BufferState& destination,
    const BufferState& source)
{
    destination.Identity = source.Identity;
    destination.Materialized = source.Materialized;
}

static LEO_FORCE_INLINE ffe_t OnePlusMultiplierLog(ffe_t skew)
{
    return LogLUT[static_cast<ffe_t>(ExpLUT[skew] ^ 1)];
}

static bool TryFourBufferButterfly(
    uint64_t buffer_bytes,
    void* a,
    void* b,
    void* c,
    void* d,
    bool inverse,
    ffe_t skew01,
    ffe_t skew23,
    ffe_t skew02)
{
    // Disabled is the public/default experiment state.  Check it before the
    // comparatively heavier kernel-resolution path so the established
    // two-way backend does not pay AVX-512 feature/translation-unit probes for
    // every complete radix-4 unit.
    const FourBufferMode four_mode =
        SelectedFourBufferMode.load(std::memory_order_relaxed);
    if (four_mode == FourBufferMode::Disabled)
        return false;

    // ResolveKernelMode performs the exact CPUID/XCR0 gate before this
    // baseline object can enter the separately targeted AVX-512 object.
    if (ResolveKernelMode() != KernelMode::Avx512Zmm ||
        !avx512four::KernelsBuilt())
    {
        return false;
    }

    const int tuple_index = avx512four::FindTupleIndex(
        skew01, skew23, skew02);
    if (tuple_index < 0 ||
        !avx512four::Apply512(
            static_cast<unsigned>(tuple_index),
            a, b, c, d, buffer_bytes, inverse,
            four_mode == FourBufferMode::Xor3))
    {
        return false;
    }

    ++LastFourBufferStatistics.FusedUnits;
    const uint64_t original_factor =
        (skew01 == kModulus ? 3 : 4) +
        (skew23 == kModulus ? 3 : 4) +
        (skew02 == kModulus ? 6 : 8);
    // Conservatively charge one fused pass as reading and writing all four
    // buffers: eight scheduled payload bytes per logical shard byte.  The
    // sentinel-bearing generated map can leave one buffer unchanged and the
    // compiler may elide those stores, so observed traffic can be lower.
    LastFourBufferStatistics.EstimatedPayloadBytesElided +=
        buffer_bytes * (original_factor - 8);
    return true;
}

bool FourBufferButterflyBufferForTesting(
    uint64_t buffer_bytes,
    void* a,
    void* b,
    void* c,
    void* d,
    bool inverse,
    ffe_t skew01,
    ffe_t skew23,
    ffe_t skew02)
{
    return TryFourBufferButterfly(
        buffer_bytes, a, b, c, d, inverse,
        skew01, skew23, skew02);
}

static void FormalDerivativeLeftUntracked(
    uint64_t buffer_bytes,
    unsigned half_count,
    void** work)
{
    for (unsigned index = 1; index < half_count; ++index)
    {
        const unsigned width = ((index ^ (index - 1)) + 1) >> 1;
        XorBuffers(
            buffer_bytes,
            width,
            work + index - width,
            work + index);
    }
}

void FormalDerivativeTopFFTForTesting(
    uint64_t buffer_bytes,
    unsigned count,
    void** work)
{
    const unsigned half_count = count >> 1;
    FormalDerivativeLeftUntracked(buffer_bytes, half_count, work);
    ApplyFormalDerivativeTopFFTBoundary(buffer_bytes, half_count, work);
}

static void IFFT_DIT_Untracked(
    uint64_t buffer_bytes,
    unsigned count_truncated,
    void** work,
    unsigned count,
    unsigned skew_base)
{
    unsigned distance = 1;
    for (; distance <= count >> 2; distance <<= 2)
    {
        const unsigned span = distance << 2;
        for (unsigned range = 0;
            range < count_truncated;
            range += span)
        {
            const ffe_t skew01 =
                FFTSkewPadded[skew_base + range + distance];
            const ffe_t skew23 =
                FFTSkewPadded[skew_base + range + distance * 3];
            const ffe_t skew02 =
                FFTSkewPadded[skew_base + range + distance * 2];
            const bool complete =
                range + distance * 2 < count_truncated;
            for (unsigned index = range; index < range + distance; ++index)
            {
                if (complete && TryFourBufferButterfly(
                        buffer_bytes,
                        work[index],
                        work[index + distance],
                        work[index + distance * 2],
                        work[index + distance * 3],
                        true,
                        skew01,
                        skew23,
                        skew02))
                {
                    continue;
                }

                // Preserve the exact two-layer IFFT ordering for every unit
                // that is truncated, narrow, unavailable, or not generated.
                IFFTButterflyBuffer(buffer_bytes,
                    work[index], work[index + distance], skew01);
                if (complete)
                {
                    IFFTButterflyBuffer(buffer_bytes,
                        work[index + distance * 2],
                        work[index + distance * 3], skew23);
                }
                IFFTButterflyBuffer(buffer_bytes,
                    work[index], work[index + distance * 2], skew02);
                IFFTButterflyBuffer(buffer_bytes,
                    work[index + distance],
                    work[index + distance * 3], skew02);
            }
        }
    }

    // A power-of-two transform has at most one unpaired layer.
    if (distance < count)
    {
        const unsigned span = distance << 1;
        for (unsigned range = 0;
            range < count_truncated;
            range += span)
        {
            const ffe_t skew = FFTSkewPadded[
                skew_base + range + distance];
            for (unsigned index = range;
                index < range + distance;
                ++index)
            {
                IFFTButterflyBuffer(buffer_bytes,
                    work[index], work[index + distance], skew);
            }
        }
    }
}

static void FFT_DIT_UntrackedFrom(
    uint64_t buffer_bytes,
    void** work,
    unsigned count_truncated,
    unsigned count,
    unsigned skew_base,
    unsigned initial_high_distance)
{
    (void)count;
    unsigned high_distance = initial_high_distance;
    for (; high_distance >= 2; high_distance >>= 2)
    {
        const unsigned distance = high_distance >> 1;
        const unsigned span = high_distance << 1;
        for (unsigned range = 0;
            range < count_truncated;
            range += span)
        {
            const ffe_t skew01 =
                FFTSkewPadded[skew_base + range + distance];
            const ffe_t skew23 =
                FFTSkewPadded[skew_base + range + distance * 3];
            const ffe_t skew02 =
                FFTSkewPadded[skew_base + range + high_distance];
            const bool complete =
                range + high_distance < count_truncated;
            for (unsigned index = range; index < range + distance; ++index)
            {
                if (complete && TryFourBufferButterfly(
                        buffer_bytes,
                        work[index],
                        work[index + distance],
                        work[index + high_distance],
                        work[index + high_distance + distance],
                        false,
                        skew01,
                        skew23,
                        skew02))
                {
                    continue;
                }

                FFTButterflyBuffer(buffer_bytes,
                    work[index], work[index + high_distance], skew02);
                FFTButterflyBuffer(buffer_bytes,
                    work[index + distance],
                    work[index + high_distance + distance], skew02);
                FFTButterflyBuffer(buffer_bytes,
                    work[index], work[index + distance], skew01);
                if (complete)
                {
                    FFTButterflyBuffer(buffer_bytes,
                        work[index + high_distance],
                        work[index + high_distance + distance], skew23);
                }
            }
        }
    }


    if (high_distance != 0)
    {
        const unsigned span = high_distance << 1;
        for (unsigned range = 0;
            range < count_truncated;
            range += span)
        {
            const ffe_t skew = FFTSkewPadded[
                skew_base + range + high_distance];
            for (unsigned index = range;
                index < range + high_distance;
                ++index)
            {
                FFTButterflyBuffer(buffer_bytes,
                    work[index], work[index + high_distance], skew);
            }
        }
    }
}

static void FFT_DIT_Untracked(
    uint64_t buffer_bytes,
    void** work,
    unsigned count_truncated,
    unsigned count,
    unsigned skew_base)
{
    FFT_DIT_UntrackedFrom(
        buffer_bytes,
        work,
        count_truncated,
        count,
        skew_base,
        count >> 1);
}

static void TrackedButterfly(
    MaterializationContext& context,
    bool inverse,
    void* x,
    void* y,
    BufferState& x_state,
    BufferState& y_state,
    ffe_t skew)
{
    const bool x_zero = IsZero(x_state);
    const bool y_zero = IsZero(y_state);

    if (x_zero && y_zero)
    {
        context.RecordSkippedButterfly(skew, false);
        return;
    }

    // Both FFT directions reduce to y ^= x for the sentinel coefficient.
    if (skew == kModulus)
    {
        if (x_zero)
        {
            context.RecordSkippedButterfly(skew, false);
            return;
        }
        if (!y_zero && x_state.Identity == y_state.Identity)
        {
            context.SetZero(y_state, false);
            context.RecordSkippedButterfly(skew, true);
            return;
        }
        if (y_zero)
        {
            memmove(y, x, static_cast<size_t>(context.BufferBytes));
            CopyBufferState(y_state, x_state);
            context.RecordReducedButterfly(
                skew, context.BufferBytes * 2);
            return;
        }

        XorBuffer(context.BufferBytes, y, x);
        context.SetFresh(y_state);
        return;
    }

    // IFFT(y == x) gives y' = 0 and x' = x for every coefficient.
    if (inverse && !x_zero &&
        x_state.Identity == y_state.Identity)
    {
        context.SetZero(y_state, false);
        context.RecordSkippedButterfly(skew, true);
        return;
    }

    // For coefficient one, FFT(x == y) gives x' = 0 and y' = y.
    if (!inverse && skew == 0 && !x_zero &&
        x_state.Identity == y_state.Identity)
    {
        context.SetZero(x_state, false);
        context.RecordSkippedButterfly(skew, true);
        return;
    }

    if (!inverse && y_zero)
    {
        // FFT(x, 0) = (x, x).
        memmove(y, x, static_cast<size_t>(context.BufferBytes));
        CopyBufferState(y_state, x_state);
        context.RecordReducedButterfly(
            skew, context.BufferBytes * 2);
        return;
    }

    if (!inverse && x_zero && skew == 0)
    {
        // FFT(0, y) with multiplier one is (y, 0).
        memmove(x, y, static_cast<size_t>(context.BufferBytes));
        CopyBufferState(x_state, y_state);
        context.SetZero(y_state, false);
        context.RecordReducedButterfly(
            skew, context.BufferBytes * 2);
        return;
    }

    if (inverse && x_zero)
    {
        // IFFT(0, y) = (M*y, y).
        MultiplyBuffer(context.BufferBytes, x, y, skew);
        if (skew == 0)
            CopyBufferState(x_state, y_state);
        else
            context.SetFresh(x_state);
        context.RecordReducedButterfly(
            skew, context.BufferBytes * 2);
        return;
    }

    if (inverse && y_zero)
    {
        // IFFT(x, 0) = ((1+M)*x, x).  Preserve x for y before applying
        // the combined coefficient in place.
        memmove(y, x, static_cast<size_t>(context.BufferBytes));
        CopyBufferState(y_state, x_state);
        const ffe_t combined_value = static_cast<ffe_t>(ExpLUT[skew] ^ 1);
        uint64_t replacement_bytes = context.BufferBytes * 2;
        if (combined_value == 0)
            context.SetZero(x_state, false);
        else
        {
            MultiplyBuffer(
                context.BufferBytes, x, x, OnePlusMultiplierLog(skew));
            context.SetFresh(x_state);
            replacement_bytes += context.BufferBytes * 2;
        }
        context.RecordReducedButterfly(skew, replacement_bytes);
        return;
    }

    // The remaining x-zero forward case needs both outputs and has no
    // no-temporary reduction using the existing kernels.  Materialize only at
    // this consumption boundary, then preserve the established butterfly.
    context.MaterializeZero(x, x_state);
    context.MaterializeZero(y, y_state);
    if (inverse)
        IFFTButterflyBuffer(context.BufferBytes, x, y, skew);
    else
        FFTButterflyBuffer(context.BufferBytes, x, y, skew);
    context.SetFresh(x_state);
    context.SetFresh(y_state);
}

void TrackedButterflyBufferForTesting(
    uint64_t buffer_bytes,
    void* x,
    void* y,
    bool inverse,
    bool x_is_zero,
    bool y_is_zero,
    bool equal_nonzero,
    ffe_t skew)
{
    MaterializationContext context(buffer_bytes);
    BufferState x_state;
    BufferState y_state;
    if (x_is_zero)
        context.SetZero(x_state, true);
    else
        context.SetFresh(x_state);
    if (y_is_zero)
        context.SetZero(y_state, true);
    else if (equal_nonzero)
        CopyBufferState(y_state, x_state);
    else
        context.SetFresh(y_state);

    TrackedButterfly(
        context, inverse, x, y, x_state, y_state, skew);
    context.MaterializeZero(x, x_state);
    context.MaterializeZero(y, y_state);
    LastMaterializationStatistics = context.Statistics;
}

static LEO_FORCE_INLINE bool FourStatesAreFusionSafe(
    const BufferState& a,
    const BufferState& b,
    const BufferState& c,
    const BufferState& d)
{
    return !IsZero(a) && !IsZero(b) && !IsZero(c) && !IsZero(d) &&
        a.Materialized && b.Materialized && c.Materialized && d.Materialized &&
        a.Identity != b.Identity && a.Identity != c.Identity &&
        a.Identity != d.Identity && b.Identity != c.Identity &&
        b.Identity != d.Identity && c.Identity != d.Identity;
}

static LEO_FORCE_INLINE void AdvanceGenericButterflyStates(
    MaterializationContext& context,
    BufferState& x,
    BufferState& y,
    ffe_t skew)
{
    // The sentinel preserves x and writes only y.  Every non-sentinel generic
    // butterfly writes two values whose equality is not statically known.
    if (skew != kModulus)
        context.SetFresh(x);
    context.SetFresh(y);
}

static void AdvanceFourBufferStates(
    MaterializationContext& context,
    bool inverse,
    BufferState& a,
    BufferState& b,
    BufferState& c,
    BufferState& d,
    ffe_t skew01,
    ffe_t skew23,
    ffe_t skew02)
{
    if (inverse)
    {
        AdvanceGenericButterflyStates(context, a, b, skew01);
        AdvanceGenericButterflyStates(context, c, d, skew23);
        AdvanceGenericButterflyStates(context, a, c, skew02);
        AdvanceGenericButterflyStates(context, b, d, skew02);
    }
    else
    {
        AdvanceGenericButterflyStates(context, a, c, skew02);
        AdvanceGenericButterflyStates(context, b, d, skew02);
        AdvanceGenericButterflyStates(context, a, b, skew01);
        AdvanceGenericButterflyStates(context, c, d, skew23);
    }
}

static void IFFT_DIT(
    uint64_t buffer_bytes,
    unsigned count_truncated,
    void** work,
    BufferState* states,
    MaterializationContext& context,
    unsigned count,
    unsigned skew_base)
{
    (void)buffer_bytes;
    unsigned distance = 1;
    for (; distance <= count >> 2; distance <<= 2)
    {
        const unsigned span = distance << 2;
        for (unsigned range = 0;
            range < count_truncated;
            range += span)
        {
            const ffe_t skew01 =
                FFTSkewPadded[skew_base + range + distance];
            const ffe_t skew23 =
                FFTSkewPadded[skew_base + range + distance * 3];
            const ffe_t skew02 =
                FFTSkewPadded[skew_base + range + distance * 2];
            const bool complete =
                range + distance * 2 < count_truncated;
            for (unsigned index = range; index < range + distance; ++index)
            {
                BufferState& a_state = states[index];
                BufferState& b_state = states[index + distance];
                BufferState& c_state = states[index + distance * 2];
                BufferState& d_state = states[index + distance * 3];
                if (complete && FourStatesAreFusionSafe(
                        a_state, b_state, c_state, d_state) &&
                    TryFourBufferButterfly(
                        context.BufferBytes,
                        work[index],
                        work[index + distance],
                        work[index + distance * 2],
                        work[index + distance * 3],
                        true,
                        skew01,
                        skew23,
                        skew02))
                {
                    AdvanceFourBufferStates(
                        context, true,
                        a_state, b_state, c_state, d_state,
                        skew01, skew23, skew02);
                    continue;
                }

                TrackedButterfly(context, true,
                    work[index], work[index + distance],
                    a_state, b_state, skew01);
                if (complete)
                {
                    TrackedButterfly(context, true,
                        work[index + distance * 2],
                        work[index + distance * 3],
                        c_state, d_state, skew23);
                }
                TrackedButterfly(context, true,
                    work[index], work[index + distance * 2],
                    a_state, c_state, skew02);
                TrackedButterfly(context, true,
                    work[index + distance],
                    work[index + distance * 3],
                    b_state, d_state, skew02);
            }
        }
    }

    if (distance < count)
    {
        const unsigned span = distance << 1;
        for (unsigned range = 0;
            range < count_truncated;
            range += span)
        {
            const ffe_t skew = FFTSkewPadded[
                skew_base + range + distance];
            for (unsigned index = range;
                index < range + distance;
                ++index)
            {
                TrackedButterfly(context, true,
                    work[index], work[index + distance],
                    states[index], states[index + distance], skew);
            }
        }
    }
}

static void FFT_DIT_From(
    uint64_t buffer_bytes,
    void** work,
    BufferState* states,
    MaterializationContext& context,
    unsigned count_truncated,
    unsigned count,
    unsigned skew_base,
    unsigned initial_high_distance)
{
    (void)buffer_bytes;
    (void)count;
    unsigned high_distance = initial_high_distance;
    for (; high_distance >= 2; high_distance >>= 2)
    {
        const unsigned distance = high_distance >> 1;
        const unsigned span = high_distance << 1;
        for (unsigned range = 0;
            range < count_truncated;
            range += span)
        {
            const ffe_t skew01 =
                FFTSkewPadded[skew_base + range + distance];
            const ffe_t skew23 =
                FFTSkewPadded[skew_base + range + distance * 3];
            const ffe_t skew02 =
                FFTSkewPadded[skew_base + range + high_distance];
            const bool complete =
                range + high_distance < count_truncated;
            for (unsigned index = range; index < range + distance; ++index)
            {
                BufferState& a_state = states[index];
                BufferState& b_state = states[index + distance];
                BufferState& c_state = states[index + high_distance];
                BufferState& d_state =
                    states[index + high_distance + distance];
                if (complete && FourStatesAreFusionSafe(
                        a_state, b_state, c_state, d_state) &&
                    TryFourBufferButterfly(
                        context.BufferBytes,
                        work[index],
                        work[index + distance],
                        work[index + high_distance],
                        work[index + high_distance + distance],
                        false,
                        skew01,
                        skew23,
                        skew02))
                {
                    AdvanceFourBufferStates(
                        context, false,
                        a_state, b_state, c_state, d_state,
                        skew01, skew23, skew02);
                    continue;
                }

                TrackedButterfly(context, false,
                    work[index], work[index + high_distance],
                    a_state, c_state, skew02);
                TrackedButterfly(context, false,
                    work[index + distance],
                    work[index + high_distance + distance],
                    b_state, d_state, skew02);
                TrackedButterfly(context, false,
                    work[index], work[index + distance],
                    a_state, b_state, skew01);
                if (complete)
                {
                    TrackedButterfly(context, false,
                        work[index + high_distance],
                        work[index + high_distance + distance],
                        c_state, d_state, skew23);
                }
            }
        }
    }

    if (high_distance != 0)
    {
        const unsigned span = high_distance << 1;
        for (unsigned range = 0;
            range < count_truncated;
            range += span)
        {
            const ffe_t skew = FFTSkewPadded[
                skew_base + range + high_distance];
            for (unsigned index = range;
                index < range + high_distance;
                ++index)
            {
                TrackedButterfly(context, false,
                    work[index], work[index + high_distance],
                    states[index], states[index + high_distance], skew);
            }
        }
    }
}

static void FFT_DIT(
    uint64_t buffer_bytes,
    void** work,
    BufferState* states,
    MaterializationContext& context,
    unsigned count_truncated,
    unsigned count,
    unsigned skew_base)
{
    FFT_DIT_From(
        buffer_bytes,
        work,
        states,
        context,
        count_truncated,
        count,
        skew_base,
        count >> 1);
}

static void TrackedXorBuffers(
    MaterializationContext& context,
    unsigned count,
    void** destination,
    BufferState* destination_states,
    void** source,
    const BufferState* source_states)
{
    void* generic_destination[kOrder];
    void* generic_source[kOrder];
    BufferState* generic_states[kOrder];
    unsigned generic_count = 0;

    for (unsigned index = 0; index < count; ++index)
    {
        BufferState& destination_state = destination_states[index];
        const BufferState& source_state = source_states[index];
        if (destination[index] == source[index])
        {
            context.SetZero(destination_state, false);
            context.RecordSkippedXor(true);
            continue;
        }
        if (IsZero(source_state))
        {
            context.RecordSkippedXor(false);
            continue;
        }
        if (IsZero(destination_state))
        {
            memmove(
                destination[index], source[index],
                static_cast<size_t>(context.BufferBytes));
            CopyBufferState(destination_state, source_state);
            ++context.Statistics.XorsReplacedByCopies;
            context.Statistics.EstimatedPayloadBytesElided +=
                static_cast<int64_t>(context.BufferBytes);
            continue;
        }
        if (destination_state.Identity == source_state.Identity)
        {
            context.SetZero(destination_state, false);
            context.RecordSkippedXor(true);
            continue;
        }

        generic_destination[generic_count] = destination[index];
        generic_source[generic_count] = source[index];
        generic_states[generic_count] = &destination_state;
        ++generic_count;
    }

    if (generic_count == 0)
        return;
    XorBuffers(
        context.BufferBytes,
        generic_count,
        generic_destination,
        generic_source);
    for (unsigned index = 0; index < generic_count; ++index)
        context.SetFresh(*generic_states[index]);
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

static void FFT_DIT_ErrorBits_UntrackedFrom(
    uint64_t buffer_bytes,
    void** work,
    unsigned count_truncated,
    unsigned count,
    unsigned skew_base,
    const ErrorBitfield& error_bits,
    unsigned initial_high_distance)
{
    (void)count;
    unsigned mip_level = initial_high_distance == 0
        ? 0
        : LastNonzeroBit32(initial_high_distance) + 1;
    unsigned high_distance = initial_high_distance;
    for (; high_distance >= 2;
        high_distance >>= 2, mip_level -= 2)
    {
        const unsigned distance = high_distance >> 1;
        const unsigned span = high_distance << 1;
        for (unsigned range = 0;
            range < count_truncated;
            range += span)
        {
            const bool complete =
                range + high_distance < count_truncated;
            const bool high_needed =
                error_bits.IsNeeded(mip_level, range);
            const bool low01_needed =
                error_bits.IsNeeded(mip_level - 1, range);
            const bool low23_needed = complete &&
                error_bits.IsNeeded(mip_level - 1,
                    range + high_distance);
            const ffe_t skew01 =
                FFTSkewPadded[skew_base + range + distance];
            const ffe_t skew23 =
                FFTSkewPadded[skew_base + range + distance * 3];
            const ffe_t skew02 =
                FFTSkewPadded[skew_base + range + high_distance];
            for (unsigned index = range; index < range + distance; ++index)
            {
                if (complete && high_needed && low01_needed && low23_needed &&
                    TryFourBufferButterfly(
                        buffer_bytes,
                        work[index],
                        work[index + distance],
                        work[index + high_distance],
                        work[index + high_distance + distance],
                        false,
                        skew01,
                        skew23,
                        skew02))
                {
                    continue;
                }

                if (high_needed)
                {
                    FFTButterflyBuffer(buffer_bytes,
                        work[index], work[index + high_distance], skew02);
                    FFTButterflyBuffer(buffer_bytes,
                        work[index + distance],
                        work[index + high_distance + distance], skew02);
                }
                if (low01_needed)
                {
                    FFTButterflyBuffer(buffer_bytes,
                        work[index], work[index + distance], skew01);
                }
                if (low23_needed)
                {
                    FFTButterflyBuffer(buffer_bytes,
                        work[index + high_distance],
                        work[index + high_distance + distance], skew23);
                }
            }
        }
    }

    if (high_distance != 0)
    {
        const unsigned span = high_distance << 1;
        for (unsigned range = 0;
            range < count_truncated;
            range += span)
        {
            if (!error_bits.IsNeeded(mip_level, range))
                continue;
            const ffe_t skew = FFTSkewPadded[
                skew_base + range + high_distance];
            for (unsigned index = range;
                index < range + high_distance;
                ++index)
            {
                FFTButterflyBuffer(buffer_bytes,
                    work[index], work[index + high_distance], skew);
            }
        }
    }
}

static void FFT_DIT_ErrorBits_From(
    uint64_t buffer_bytes,
    void** work,
    BufferState* states,
    MaterializationContext& context,
    unsigned count_truncated,
    unsigned count,
    unsigned skew_base,
    const ErrorBitfield& error_bits,
    unsigned initial_high_distance)
{
    (void)buffer_bytes;
    (void)count;
    unsigned mip_level = initial_high_distance == 0
        ? 0
        : LastNonzeroBit32(initial_high_distance) + 1;
    unsigned high_distance = initial_high_distance;
    for (; high_distance >= 2;
        high_distance >>= 2, mip_level -= 2)
    {
        const unsigned distance = high_distance >> 1;
        const unsigned span = high_distance << 1;
        for (unsigned range = 0;
            range < count_truncated;
            range += span)
        {
            const bool complete =
                range + high_distance < count_truncated;
            const bool high_needed =
                error_bits.IsNeeded(mip_level, range);
            const bool low01_needed =
                error_bits.IsNeeded(mip_level - 1, range);
            const bool low23_needed = complete &&
                error_bits.IsNeeded(mip_level - 1,
                    range + high_distance);
            const ffe_t skew01 =
                FFTSkewPadded[skew_base + range + distance];
            const ffe_t skew23 =
                FFTSkewPadded[skew_base + range + distance * 3];
            const ffe_t skew02 =
                FFTSkewPadded[skew_base + range + high_distance];
            for (unsigned index = range; index < range + distance; ++index)
            {
                BufferState& a_state = states[index];
                BufferState& b_state = states[index + distance];
                BufferState& c_state = states[index + high_distance];
                BufferState& d_state =
                    states[index + high_distance + distance];
                if (complete && high_needed && low01_needed && low23_needed &&
                    FourStatesAreFusionSafe(
                        a_state, b_state, c_state, d_state) &&
                    TryFourBufferButterfly(
                        context.BufferBytes,
                        work[index],
                        work[index + distance],
                        work[index + high_distance],
                        work[index + high_distance + distance],
                        false,
                        skew01,
                        skew23,
                        skew02))
                {
                    AdvanceFourBufferStates(
                        context, false,
                        a_state, b_state, c_state, d_state,
                        skew01, skew23, skew02);
                    continue;
                }

                if (high_needed)
                {
                    TrackedButterfly(context, false,
                        work[index], work[index + high_distance],
                        a_state, c_state, skew02);
                    TrackedButterfly(context, false,
                        work[index + distance],
                        work[index + high_distance + distance],
                        b_state, d_state, skew02);
                }
                if (low01_needed)
                {
                    TrackedButterfly(context, false,
                        work[index], work[index + distance],
                        a_state, b_state, skew01);
                }
                if (low23_needed)
                {
                    TrackedButterfly(context, false,
                        work[index + high_distance],
                        work[index + high_distance + distance],
                        c_state, d_state, skew23);
                }
            }
        }
    }

    if (high_distance != 0)
    {
        const unsigned span = high_distance << 1;
        for (unsigned range = 0;
            range < count_truncated;
            range += span)
        {
            if (!error_bits.IsNeeded(mip_level, range))
                continue;
            const ffe_t skew = FFTSkewPadded[
                skew_base + range + high_distance];
            for (unsigned index = range;
                index < range + high_distance;
                ++index)
            {
                TrackedButterfly(context, false,
                    work[index], work[index + high_distance],
                    states[index], states[index + high_distance], skew);
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

struct LocatorShiftTerm
{
    uint8_t BaseLog;
    uint8_t Inverse;
};

// GCC's -O3 -floop-unroll-and-jam pass combines adjacent terms in the nested
// accumulator below and turns the 256-shift loop into scalar stack updates.
// Keeping that single pass off lets the existing vectorizer widen each term's
// independent uint8_t-to-uint16_t accumulation.  Other compilers either do not
// expose this pass through a function attribute or already generate the
// intended vector loop.  Apply the same scoped guard to native and MinGW GCC;
// MSVC and Clang do not enter this compiler-specific path.
#if defined(__GNUC__) && !defined(__clang__) && (__GNUC__ >= 8)
    #define LEO_FF8XOR_NO_UNROLL_AND_JAM \
        __attribute__((optimize("no-loop-unroll-and-jam")))
#else
    #define LEO_FF8XOR_NO_UNROLL_AND_JAM
#endif

static_assert(
    kOrder * generated::kMultiplyMaxGateCount <= UINT16_MAX,
    "locator gate totals fit in uint16_t");
static_assert(
    kOrder * generated::kMultiplyMaxDepth <= UINT16_MAX,
    "locator depth totals fit in uint16_t");

static LEO_FF8XOR_NO_UNROLL_AND_JAM unsigned SelectLocatorShiftTerms(
    const LocatorShiftTerm* terms,
    unsigned term_count)
{
    LEO_ALIGNED uint16_t gate_totals[kOrder] = {};

    // Accumulate all 255 candidate rotations term by term.  The generated
    // rows contain one extra duplicate entry so this hot inner loop has the
    // power-of-two trip count 256; gate_totals[255] is intentionally ignored.
    for (unsigned term_index = 0; term_index < term_count; ++term_index)
    {
        const LocatorShiftTerm& term = terms[term_index];
        const uint8_t* row = term.Inverse
            ? generated::kLocatorInverseGateCosts
            : generated::kLocatorPositiveGateCosts;
        row += term.BaseLog;
        for (unsigned shift = 0; shift < kOrder; ++shift)
        {
            gate_totals[shift] = static_cast<uint16_t>(
                gate_totals[shift] + row[shift]);
        }
    }

    uint16_t best_gates = UINT16_MAX;
    for (unsigned shift = 0; shift < kModulus; ++shift)
    {
        if (gate_totals[shift] < best_gates)
            best_gates = gate_totals[shift];
    }

    unsigned best_shift = 0;
    unsigned best_depth = UINT_MAX;
    for (unsigned shift = 0; shift < kModulus; ++shift)
    {
        if (gate_totals[shift] != best_gates)
            continue;

        unsigned depth = 0;
        for (unsigned term_index = 0; term_index < term_count; ++term_index)
        {
            const LocatorShiftTerm& term = terms[term_index];
            const uint8_t* row = term.Inverse
                ? generated::kLocatorInverseDepthCosts
                : generated::kLocatorPositiveDepthCosts;
            depth += row[term.BaseLog + shift];
        }

        // Strict improvement plus ascending scan preserves the old lowest
        // numeric shift tie-break when both gate count and depth are equal.
        if (depth < best_depth)
        {
            best_depth = depth;
            best_shift = shift;
        }
    }

    return best_shift;
}

#undef LEO_FF8XOR_NO_UNROLL_AND_JAM

#ifdef _MSC_VER
    #define LEO_FF8XOR_LOCATOR_NOINLINE __declspec(noinline)
#elif defined(__GNUC__) || defined(__clang__)
    #define LEO_FF8XOR_LOCATOR_NOINLINE __attribute__((noinline))
#else
    #define LEO_FF8XOR_LOCATOR_NOINLINE
#endif

LEO_FF8XOR_LOCATOR_NOINLINE unsigned SelectLocatorShiftForTesting(
    const ffe_t* logarithms,
    const bool* inverse,
    unsigned count)
{
    if (count > kOrder ||
        (count != 0 && (logarithms == NULL || inverse == NULL)))
    {
        return kModulus;
    }

    LocatorShiftTerm terms[kOrder];
    for (unsigned index = 0; index < count; ++index)
    {
        terms[index].BaseLog = static_cast<uint8_t>(
            CanonicalLog(logarithms[index]));
        terms[index].Inverse = inverse[index] ? 1 : 0;
    }
    return SelectLocatorShiftTerms(terms, count);
}

LEO_FF8XOR_LOCATOR_NOINLINE unsigned
SelectLocatorShiftReferenceForTesting(
    const ffe_t* logarithms,
    const bool* inverse,
    unsigned count)
{
    if (count > kOrder ||
        (count != 0 && (logarithms == NULL || inverse == NULL)))
    {
        return kModulus;
    }

    uint64_t best_gates = UINT64_MAX;
    uint64_t best_depth = UINT64_MAX;
    unsigned best_shift = 0;
    for (unsigned shift = 0; shift < kModulus; ++shift)
    {
        uint64_t gates = 0;
        uint64_t depth = 0;
        for (unsigned index = 0; index < count; ++index)
        {
            const ffe_t log = inverse[index]
                ? InverseShiftedLog(logarithms[index], shift)
                : ShiftedLog(logarithms[index], shift);
            gates += generated::kMultiplyGateCounts[log];
            depth += generated::kMultiplyDepths[log];
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

#undef LEO_FF8XOR_LOCATOR_NOINLINE

static unsigned SelectLocatorShift(
    const ffe_t* error_locations,
    unsigned original_count,
    unsigned recovery_count,
    unsigned m,
    const void* const* original,
    const void* const* recovery)
{
    LocatorShiftTerm terms[kOrder];
    unsigned term_count = 0;

    for (unsigned index = 0; index < recovery_count; ++index)
    {
        if (recovery[index])
        {
            terms[term_count].BaseLog = static_cast<uint8_t>(
                CanonicalLog(error_locations[index]));
            terms[term_count].Inverse = 0;
            ++term_count;
        }
    }

    for (unsigned index = 0; index < original_count; ++index)
    {
        terms[term_count].BaseLog = static_cast<uint8_t>(
            CanonicalLog(error_locations[m + index]));
        terms[term_count].Inverse = original[index] ? 0 : 1;
        ++term_count;
    }

    return SelectLocatorShiftTerms(terms, term_count);
}


//------------------------------------------------------------------------------
// Reed-Solomon encoding

static const uint64_t kMaterializationMinimumBufferBytes = 65536;

static void ReedSolomonEncodeUntracked(
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

        IFFT_DIT_Untracked(
            buffer_bytes,
            chunk_count,
            chunk_work,
            m,
            m + chunk_start);
        if (chunk_start != 0)
            XorBuffers(buffer_bytes, m, work, chunk_work);
        chunk_start += m;
    }
    FFT_DIT_Untracked(buffer_bytes, work, recovery_count, m, 0);
}

void ReedSolomonEncode(
    uint64_t buffer_bytes,
    unsigned original_count,
    unsigned recovery_count,
    unsigned m,
    const void* const* data,
    void** work)
{
    ResetFourBufferStatistics();
    // Full chunks contain no structural zeros or copy identities, so tracking
    // can only add schedule overhead.  Small payloads likewise cannot amortize
    // the state machine; preserve the established schedule in both cases.
    if (buffer_bytes < kMaterializationMinimumBufferBytes ||
        original_count % m == 0)
    {
        ReedSolomonEncodeUntracked(
            buffer_bytes,
            original_count,
            recovery_count,
            m,
            data,
            work);
        const MaterializationStatistics empty = {};
        LastMaterializationStatistics = empty;
        return;
    }

    MaterializationContext context(buffer_bytes);
    BufferState accumulated_states[kOrder];
    BufferState temporary_states[kOrder];
    unsigned chunk_start = 0;
    while (chunk_start < original_count)
    {
        const unsigned remaining = original_count - chunk_start;
        const unsigned chunk_count = remaining < m ? remaining : m;
        void** chunk_work = chunk_start == 0 ? work : work + m;
        BufferState* chunk_states = chunk_start == 0
            ? accumulated_states
            : temporary_states;

        for (unsigned index = 0; index < chunk_count; ++index)
        {
            memcpy(
                chunk_work[index],
                data[chunk_start + index],
                static_cast<size_t>(buffer_bytes));
            context.SetFresh(chunk_states[index]);
        }
        for (unsigned index = chunk_count; index < m; ++index)
            context.DeferZero(chunk_states[index]);

        // Padded index corresponding to logical FFTSkew + m - 1 + chunk_start.
        IFFT_DIT(
            buffer_bytes,
            chunk_count,
            chunk_work,
            chunk_states,
            context,
            m,
            m + chunk_start);

        if (chunk_start != 0)
        {
            TrackedXorBuffers(
                context,
                m,
                work,
                accumulated_states,
                chunk_work,
                chunk_states);
        }

        chunk_start += m;
    }

    FFT_DIT(
        buffer_bytes,
        work,
        accumulated_states,
        context,
        recovery_count,
        m,
        0);
    for (unsigned index = 0; index < recovery_count; ++index)
        context.MaterializeZero(work[index], accumulated_states[index]);
    LastMaterializationStatistics = context.Statistics;
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
    ResetFourBufferStatistics();
    MaterializationContext context(buffer_bytes);
    BufferState states[kOrder];
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

    if (buffer_bytes < kMaterializationMinimumBufferBytes)
    {
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
                    ShiftedLog(
                        error_locations[m + index], locator_shift));
            }
            else
                memset(
                    work[m + index], 0,
                    static_cast<size_t>(buffer_bytes));
        }
        for (unsigned index = m + original_count; index < n; ++index)
            memset(work[index], 0, static_cast<size_t>(buffer_bytes));

        IFFT_DIT_Untracked(
            buffer_bytes, m + original_count, work, n, 0);
        const unsigned half_count = n >> 1;
        LEO_DEBUG_ASSERT(FFTSkewPadded[half_count] == kModulus);
        FormalDerivativeLeftUntracked(
            buffer_bytes, half_count, work);
        ApplyFormalDerivativeTopFFTBoundary(
            buffer_bytes, half_count, work);

        const unsigned output_count = m + original_count;
#ifdef LEO_ERROR_BITFIELD_OPT
        FFT_DIT_ErrorBits_UntrackedFrom(
            buffer_bytes,
            work,
            output_count,
            n,
            0,
            error_bits,
            n >> 2);
#else
        FFT_DIT_UntrackedFrom(
            buffer_bytes, work, output_count, n, 0, n >> 2);
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
        const MaterializationStatistics empty = {};
        LastMaterializationStatistics = empty;
        return;
    }

    for (unsigned index = 0; index < recovery_count; ++index)
    {
        if (recovery[index])
        {
            MultiplyBuffer(
                buffer_bytes,
                work[index],
                recovery[index],
                ShiftedLog(error_locations[index], locator_shift));
            context.SetFresh(states[index]);
        }
        else
            context.DeferZero(states[index]);
    }
    for (unsigned index = recovery_count; index < m; ++index)
        context.DeferZero(states[index]);

    for (unsigned index = 0; index < original_count; ++index)
    {
        if (original[index])
        {
            MultiplyBuffer(
                buffer_bytes,
                work[m + index],
                original[index],
                ShiftedLog(error_locations[m + index], locator_shift));
            context.SetFresh(states[m + index]);
        }
        else
            context.DeferZero(states[m + index]);
    }
    for (unsigned index = m + original_count; index < n; ++index)
        context.DeferZero(states[index]);

    IFFT_DIT(
        buffer_bytes,
        m + original_count,
        work,
        states,
        context,
        n,
        0);

    const unsigned half_count = n >> 1;
    for (unsigned index = 1; index < half_count; ++index)
    {
        const unsigned width = ((index ^ (index - 1)) + 1) >> 1;
        TrackedXorBuffers(
            context,
            width,
            work + index - width,
            states + index - width,
            work + index,
            states + index);
    }

    // The direct boundary reads every logical A/R input at least once.  Any
    // deferred zero must therefore be materialized before the static row
    // kernels enter their payload loops.  Every output is conservatively given
    // a fresh identity afterward; later ErrorBitfield pruning and zero handling
    // remain valid without assuming an unproved symbolic cancellation.
    for (unsigned index = 0; index < n; ++index)
        context.MaterializeZero(work[index], states[index]);
    LEO_DEBUG_ASSERT(FFTSkewPadded[half_count] == kModulus);
    ApplyFormalDerivativeTopFFTBoundary(
        buffer_bytes, half_count, work);
    for (unsigned index = 0; index < n; ++index)
        context.SetFresh(states[index]);

    const unsigned output_count = m + original_count;
#ifdef LEO_ERROR_BITFIELD_OPT
    FFT_DIT_ErrorBits_From(
        buffer_bytes,
        work,
        states,
        context,
        output_count,
        n,
        0,
        error_bits,
        n >> 2);
#else
    FFT_DIT_From(
        buffer_bytes,
        work,
        states,
        context,
        output_count,
        n,
        0,
        n >> 2);
#endif

    for (unsigned index = 0; index < original_count; ++index)
    {
        if (!original[index])
        {
            BufferState& destination_state = states[index];
            const BufferState& source_state = states[index + m];
            if (IsZero(source_state))
            {
                context.SetZero(destination_state, false);
                context.MaterializeZero(work[index], destination_state);
            }
            else
            {
                const ffe_t log = InverseShiftedLog(
                    error_locations[index + m], locator_shift);
                MultiplyBuffer(
                    buffer_bytes,
                    work[index],
                    work[index + m],
                    log);
                if (log == 0 || log == kModulus)
                    CopyBufferState(destination_state, source_state);
                else
                    context.SetFresh(destination_state);
            }
        }
    }
    LastMaterializationStatistics = context.Statistics;
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
