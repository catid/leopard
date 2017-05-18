/*
    S.-J. Lin,  T. Y. Al-Naffouri, Y. S. Han, and W.-H. Chung,
    "Novel Polynomial Basis with Fast Fourier Transform and Its Application to Reed-Solomon Erasure Codes"
    IEEE Trans. on Information Theory, pp. 6284-6299, November, 2016.
    http://ct.ee.ntust.edu.tw/it2016-2.pdf
*/

#include <string.h>
#include <time.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>


//------------------------------------------------------------------------------
// Debug

// Some bugs only repro in release mode, so this can be helpful
//#define LHC_DEBUG_IN_RELEASE

#if defined(_DEBUG) || defined(DEBUG) || defined(LHC_DEBUG_IN_RELEASE)
    #define LHC_DEBUG
    #ifdef _WIN32
        #define LHC_DEBUG_BREAK __debugbreak()
    #else
        #define LHC_DEBUG_BREAK __builtin_trap()
    #endif
    #define LHC_DEBUG_ASSERT(cond) { if (!(cond)) { LHC_DEBUG_BREAK; } }
#else
    #define LHC_DEBUG_BREAK ;
    #define LHC_DEBUG_ASSERT(cond) ;
#endif


//------------------------------------------------------------------------------
// Platform/Architecture

#if defined(ANDROID) || defined(IOS)
    #define LHC_TARGET_MOBILE
#endif // ANDROID

#if defined(__AVX2__) || (defined (_MSC_VER) && _MSC_VER >= 1900)
    #define LHC_TRY_AVX2 /* 256-bit */
    #include <immintrin.h>
    #define LHC_ALIGN_BYTES 32
#else // __AVX2__
    #define LHC_ALIGN_BYTES 16
#endif // __AVX2__

#if !defined(LHC_TARGET_MOBILE)
    // Note: MSVC currently only supports SSSE3 but not AVX2
    #include <tmmintrin.h> // SSSE3: _mm_shuffle_epi8
    #include <emmintrin.h> // SSE2
#endif // LHC_TARGET_MOBILE

#if defined(HAVE_ARM_NEON_H)
    #include <arm_neon.h>
#endif // HAVE_ARM_NEON_H

#if defined(LHC_TARGET_MOBILE)

    #define LHC_ALIGNED_ACCESSES /* Inputs must be aligned to LHC_ALIGN_BYTES */

# if defined(HAVE_ARM_NEON_H)
    // Compiler-specific 128-bit SIMD register keyword
    #define LHC_M128 uint8x16_t
    #define LHC_TRY_NEON
#else
    #define LHC_M128 uint64_t
# endif

#else // LHC_TARGET_MOBILE

    // Compiler-specific 128-bit SIMD register keyword
    #define LHC_M128 __m128i

#endif // LHC_TARGET_MOBILE

#ifdef LHC_TRY_AVX2
    // Compiler-specific 256-bit SIMD register keyword
    #define LHC_M256 __m256i
#endif

// Compiler-specific C++11 restrict keyword
#define LHC_RESTRICT __restrict

// Compiler-specific force inline keyword
#ifdef _MSC_VER
    #define LHC_FORCE_INLINE inline __forceinline
#else
    #define LHC_FORCE_INLINE inline __attribute__((always_inline))
#endif

// Compiler-specific alignment keyword
// Note: Alignment only matters for ARM NEON where it should be 16
#ifdef _MSC_VER
    #define LHC_ALIGNED __declspec(align(LHC_ALIGN_BYTES))
#else // _MSC_VER
    #define LHC_ALIGNED __attribute__((aligned(LHC_ALIGN_BYTES)))
#endif // _MSC_VER


//------------------------------------------------------------------------------
// Runtime CPU Architecture Check
//
// Feature checks stolen shamelessly from
// https://github.com/jedisct1/libsodium/blob/master/src/libsodium/sodium/runtime.c

#if defined(HAVE_ANDROID_GETCPUFEATURES)
    #include <cpu-features.h>
#endif

#if defined(LHC_TRY_NEON)
# if defined(IOS) && defined(__ARM_NEON__)
        // Requires iPhone 5S or newer
        static const bool CpuHasNeon = true;
        static const bool CpuHasNeon64 = true;
# else
        // Remember to add LOCAL_STATIC_LIBRARIES := cpufeatures
        static bool CpuHasNeon = false; // V6 / V7
        static bool CpuHasNeon64 = false; // 64-bit
# endif
#endif


#if !defined(LHC_TARGET_MOBILE)

#ifdef _MSC_VER
    #include <intrin.h> // __cpuid
    #pragma warning(disable: 4752) // found Intel(R) Advanced Vector Extensions; consider using /arch:AVX
#endif

#ifdef LHC_TRY_AVX2
static bool CpuHasAVX2 = false;
#endif
static bool CpuHasSSSE3 = false;

#define CPUID_EBX_AVX2    0x00000020
#define CPUID_ECX_SSSE3   0x00000200

static void _cpuid(unsigned int cpu_info[4U], const unsigned int cpu_info_type)
{
#if defined(_MSC_VER) && (defined(_M_X64) || defined(_M_AMD64) || defined(_M_IX86))
    __cpuid((int *) cpu_info, cpu_info_type);
#else //if defined(HAVE_CPUID)
    cpu_info[0] = cpu_info[1] = cpu_info[2] = cpu_info[3] = 0;
# ifdef __i386__
    __asm__ __volatile__ ("pushfl; pushfl; "
                          "popl %0; "
                          "movl %0, %1; xorl %2, %0; "
                          "pushl %0; "
                          "popfl; pushfl; popl %0; popfl" :
                          "=&r" (cpu_info[0]), "=&r" (cpu_info[1]) :
                          "i" (0x200000));
    if (((cpu_info[0] ^ cpu_info[1]) & 0x200000) == 0) {
        return; /* LCOV_EXCL_LINE */
    }
# endif
# ifdef __i386__
    __asm__ __volatile__ ("xchgl %%ebx, %k1; cpuid; xchgl %%ebx, %k1" :
                          "=a" (cpu_info[0]), "=&r" (cpu_info[1]),
                          "=c" (cpu_info[2]), "=d" (cpu_info[3]) :
                          "0" (cpu_info_type), "2" (0U));
# elif defined(__x86_64__)
    __asm__ __volatile__ ("xchgq %%rbx, %q1; cpuid; xchgq %%rbx, %q1" :
                          "=a" (cpu_info[0]), "=&r" (cpu_info[1]),
                          "=c" (cpu_info[2]), "=d" (cpu_info[3]) :
                          "0" (cpu_info_type), "2" (0U));
# else
    __asm__ __volatile__ ("cpuid" :
                          "=a" (cpu_info[0]), "=b" (cpu_info[1]),
                          "=c" (cpu_info[2]), "=d" (cpu_info[3]) :
                          "0" (cpu_info_type), "2" (0U));
# endif
#endif
}

#endif // defined(LHC_TARGET_MOBILE)


static void lhc_architecture_init()
{
#if defined(LHC_TRY_NEON) && defined(HAVE_ANDROID_GETCPUFEATURES)
    AndroidCpuFamily family = android_getCpuFamily();
    if (family == ANDROID_CPU_FAMILY_ARM)
    {
        if (android_getCpuFeatures() & ANDROID_CPU_ARM_FEATURE_NEON)
            CpuHasNeon = true;
    }
    else if (family == ANDROID_CPU_FAMILY_ARM64)
    {
        CpuHasNeon = true;
        if (android_getCpuFeatures() & ANDROID_CPU_ARM64_FEATURE_ASIMD)
            CpuHasNeon64 = true;
    }
#endif

#if !defined(LHC_TARGET_MOBILE)
    unsigned int cpu_info[4];

    _cpuid(cpu_info, 1);
    CpuHasSSSE3 = ((cpu_info[2] & CPUID_ECX_SSSE3) != 0);

#if defined(LHC_TRY_AVX2)
    _cpuid(cpu_info, 7);
    CpuHasAVX2 = ((cpu_info[1] & CPUID_EBX_AVX2) != 0);
#endif // LHC_TRY_AVX2

#endif // LHC_TARGET_MOBILE
}


//------------------------------------------------------------------------------
// SIMD-Safe Aligned Memory Allocations

static const unsigned kAlignmentBytes = LHC_ALIGN_BYTES;

LHC_FORCE_INLINE unsigned NextAlignedOffset(unsigned offset)
{
    return (offset + kAlignmentBytes - 1) & ~(kAlignmentBytes - 1);
}

static LHC_FORCE_INLINE uint8_t* SIMDSafeAllocate(size_t size)
{
    uint8_t* data = (uint8_t*)calloc(1, kAlignmentBytes + size);
    if (!data)
        return nullptr;
    unsigned offset = (unsigned)((uintptr_t)data % kAlignmentBytes);
    data += kAlignmentBytes - offset;
    data[-1] = (uint8_t)offset;
    return data;
}

static LHC_FORCE_INLINE void SIMDSafeFree(void* ptr)
{
    if (!ptr)
        return;
    uint8_t* data = (uint8_t*)ptr;
    unsigned offset = data[-1];
    if (offset >= kAlignmentBytes)
    {
        LHC_DEBUG_BREAK; // Should never happen
        return;
    }
    data -= kAlignmentBytes - offset;
    free(data);
}


//------------------------------------------------------------------------------
// Field

#if 0
typedef uint8_t GFSymbol;
static const unsigned kGFBits = 8;
static const unsigned kGFPolynomial = 0x11D;
GFSymbol kGFCantorBasis[kGFBits] = {
    1, 214, 152, 146, 86, 200, 88, 230
};
#else
typedef uint16_t GFSymbol;
static const unsigned kGFBits = 16;
static const unsigned kGFPolynomial = 0x1002D;
GFSymbol kGFCantorBasis[kGFBits] = {
    0x0001, 0xACCA, 0x3C0E, 0x163E,
    0xC582, 0xED2E, 0x914C, 0x4012,
    0x6C98, 0x10D8, 0x6A72, 0xB900,
    0xFDB8, 0xFB34, 0xFF38, 0x991E
};
#endif

/*
    Cantor Basis introduced by the paper:
    D. G. Cantor, "On arithmetical algorithms over finite fields",
    Journal of Combinatorial Theory, Series A, vol. 50, no. 2, pp. 285-300, 1989.
*/

static const unsigned kFieldSize = (unsigned)1 << kGFBits; //Field size
static const unsigned kFieldModulus = kFieldSize - 1;

static GFSymbol GFLog[kFieldSize];
static GFSymbol GFExp[kFieldSize];

// Initialize GFLog[], GFExp[]
static void InitField()
{
    // Use GFExp temporarily to store the monomial basis logarithm table
    GFSymbol* MonoLog = GFExp;
    unsigned state = 1;
    for (unsigned i = 0; i < kFieldModulus; ++i)
    {
        MonoLog[state] = static_cast<GFSymbol>(i);
        state <<= 1;
        if (state >= kFieldSize)
            state ^= kGFPolynomial;
    }
    MonoLog[0] = kFieldModulus;

    // Conversion to polynomial basis:

    GFLog[0] = 0;
    for (unsigned i = 0; i < kGFBits; ++i)
    {
        const GFSymbol basis = kGFCantorBasis[i];
        const unsigned width = (unsigned)(1UL << i);

        for (unsigned j = 0; j < width; ++j)
            GFLog[j + width] = GFLog[j] ^ basis;
    }

    for (unsigned i = 0; i < kFieldSize; ++i)
        GFLog[i] = MonoLog[GFLog[i]];

    for (unsigned i = 0; i < kFieldSize; ++i)
        GFExp[GFLog[i]] = i;

    GFExp[kFieldModulus] = GFExp[0];
}


//------------------------------------------------------------------------------
// Mod Q Field Operations
//
// Q is the maximum symbol value, e.g. 255 or 65535.

// z = x + y (mod Q)
static inline GFSymbol AddModQ(GFSymbol a, GFSymbol b)
{
    const unsigned sum = (unsigned)a + b;

    // Partial reduction step, allowing for Q to be returned
    return static_cast<GFSymbol>(sum + (sum >> kGFBits));
}

// z = x - y (mod Q)
static inline GFSymbol SubModQ(GFSymbol a, GFSymbol b)
{
    const unsigned dif = (unsigned)a - b;

    // Partial reduction step, allowing for Q to be returned
    return static_cast<GFSymbol>(dif + (dif >> kGFBits));
}

// vx[] += vy[] * z
static void muladd_mem(GFSymbol * LHC_RESTRICT vx, const GFSymbol * LHC_RESTRICT vy, GFSymbol z, unsigned symbolCount)
{
    for (unsigned i = 0; i < symbolCount; ++i)
    {
        const GFSymbol a = vy[i];
        if (a == 0)
            continue;

        /*
            This function consumes all the runtime.

            I have no idea how to speed this up because the GFLog and GFExp do not have any structure.
        */
        const GFSymbol sum = static_cast<GFSymbol>(AddModQ(GFLog[a], z));

        vx[i] ^= GFExp[sum];
    }
}

// return a*GFExp[b] over GF(2^r)
static GFSymbol mulE(GFSymbol a, GFSymbol b)
{
    if (a == 0)
        return 0;

    const GFSymbol sum = static_cast<GFSymbol>(AddModQ(GFLog[a], b));
    return GFExp[sum];
}

//------------------------------------------------------------------------------
// Fast Walsh-Hadamard Transform (FWHT) Mod Q
//
// Q is the maximum symbol value, e.g. 255 or 65535.

// Define this to enable the optimized version of FWHT()
#define LHC_FWHT_OPTIMIZED

#ifndef LHC_FWHT_OPTIMIZED

// Reference implementation
static void FWHT(GFSymbol* data, const unsigned bits)
{
    const unsigned size = (unsigned)(1UL << bits);
    for (unsigned width = 1; width < size; width <<= 1)
        for (unsigned i = 0; i < size; i += (width << 1))
            for (unsigned j = i; j < (width + i); ++j)
                CrossAddSubModQ(data[j], data[j + width]);
}

#else

// {a, b} = {a + b, a - b} (mod Q)
static inline void FWHT_2(GFSymbol& a, GFSymbol& b)
{
    const GFSymbol dif = SubModQ(a, b);
    const GFSymbol sum = AddModQ(a, b);
    a = sum, b = dif;
}

static inline void FWHT_4(GFSymbol* data)
{
    GFSymbol t0 = data[0];
    GFSymbol t1 = data[1];
    GFSymbol t2 = data[2];
    GFSymbol t3 = data[3];
    FWHT_2(t0, t1);
    FWHT_2(t2, t3);
    FWHT_2(t0, t2);
    FWHT_2(t1, t3);
    data[0] = t0;
    data[1] = t1;
    data[2] = t2;
    data[3] = t3;
}

static inline void FWHT_4(GFSymbol* data, unsigned s)
{
    unsigned x = 0;
    GFSymbol t0 = data[x];  x += s;
    GFSymbol t1 = data[x];  x += s;
    GFSymbol t2 = data[x];  x += s;
    GFSymbol t3 = data[x];
    FWHT_2(t0, t1);
    FWHT_2(t2, t3);
    FWHT_2(t0, t2);
    FWHT_2(t1, t3);
    unsigned y = 0;
    data[y] = t0;  y += s;
    data[y] = t1;  y += s;
    data[y] = t2;  y += s;
    data[y] = t3;
}

static inline void FWHT_8(GFSymbol* data)
{
    GFSymbol t0 = data[0];
    GFSymbol t1 = data[1];
    GFSymbol t2 = data[2];
    GFSymbol t3 = data[3];
    GFSymbol t4 = data[4];
    GFSymbol t5 = data[5];
    GFSymbol t6 = data[6];
    GFSymbol t7 = data[7];
    FWHT_2(t0, t1);
    FWHT_2(t2, t3);
    FWHT_2(t4, t5);
    FWHT_2(t6, t7);
    FWHT_2(t0, t2);
    FWHT_2(t1, t3);
    FWHT_2(t4, t6);
    FWHT_2(t5, t7);
    FWHT_2(t0, t4);
    FWHT_2(t1, t5);
    FWHT_2(t2, t6);
    FWHT_2(t3, t7);
    data[0] = t0;
    data[1] = t1;
    data[2] = t2;
    data[3] = t3;
    data[4] = t4;
    data[5] = t5;
    data[6] = t6;
    data[7] = t7;
}

static inline void FWHT_16(GFSymbol* data)
{
    GFSymbol t0 = data[0];
    GFSymbol t1 = data[1];
    GFSymbol t2 = data[2];
    GFSymbol t3 = data[3];
    GFSymbol t4 = data[4];
    GFSymbol t5 = data[5];
    GFSymbol t6 = data[6];
    GFSymbol t7 = data[7];
    GFSymbol t8 = data[8];
    GFSymbol t9 = data[9];
    GFSymbol t10 = data[10];
    GFSymbol t11 = data[11];
    GFSymbol t12 = data[12];
    GFSymbol t13 = data[13];
    GFSymbol t14 = data[14];
    GFSymbol t15 = data[15];
    FWHT_2(t0, t1);
    FWHT_2(t2, t3);
    FWHT_2(t4, t5);
    FWHT_2(t6, t7);
    FWHT_2(t8, t9);
    FWHT_2(t10, t11);
    FWHT_2(t12, t13);
    FWHT_2(t14, t15);
    FWHT_2(t0, t2);
    FWHT_2(t1, t3);
    FWHT_2(t4, t6);
    FWHT_2(t5, t7);
    FWHT_2(t8, t10);
    FWHT_2(t9, t11);
    FWHT_2(t12, t14);
    FWHT_2(t13, t15);
    FWHT_2(t0, t4);
    FWHT_2(t1, t5);
    FWHT_2(t2, t6);
    FWHT_2(t3, t7);
    FWHT_2(t8, t12);
    FWHT_2(t9, t13);
    FWHT_2(t10, t14);
    FWHT_2(t11, t15);
    FWHT_2(t0, t8);
    FWHT_2(t1, t9);
    FWHT_2(t2, t10);
    FWHT_2(t3, t11);
    FWHT_2(t4, t12);
    FWHT_2(t5, t13);
    FWHT_2(t6, t14);
    FWHT_2(t7, t15);
    data[0] = t0;
    data[1] = t1;
    data[2] = t2;
    data[3] = t3;
    data[4] = t4;
    data[5] = t5;
    data[6] = t6;
    data[7] = t7;
    data[8] = t8;
    data[9] = t9;
    data[10] = t10;
    data[11] = t11;
    data[12] = t12;
    data[13] = t13;
    data[14] = t14;
    data[15] = t15;
}

void FWHT_SmallData(GFSymbol* data, unsigned ldn)
{
    const unsigned n = (1UL << ldn);

    if (n <= 2)
    {
        if (n == 2)
            FWHT_2(data[0], data[1]);
        return;
    }

    for (unsigned ldm = ldn; ldm > 3; ldm -= 2)
    {
        unsigned m = (1UL << ldm);
        unsigned m4 = (m >> 2);
        for (unsigned r = 0; r < n; r += m)
            for (unsigned j = 0; j < m4; j++)
                FWHT_4(data + j + r, m4);
    }

    if (ldn & 1)
    {
        for (unsigned i0 = 0; i0 < n; i0 += 8)
            FWHT_8(data + i0);
    }
    else
    {
        for (unsigned i0 = 0; i0 < n; i0 += 4)
            FWHT_4(data + i0);
    }
}

// Decimation in time version of the transform
static void FWHT(GFSymbol* data, const unsigned ldn)
{
    if (ldn <= 13)
    {
        FWHT_SmallData(data, ldn);
        return;
    }

    FWHT_2(data[2], data[3]);
    FWHT_4(data + 4);
    FWHT_8(data + 8);
    FWHT_16(data + 16);
    for (unsigned ldm = 5; ldm < ldn; ++ldm)
        FWHT(data + (unsigned)(1UL << ldm), ldm);

    for (unsigned ldm = 1; ldm <= ldn; ++ldm)
    {
        const unsigned m = (1UL << ldm);
        const unsigned mh = (m >> 1);
        for (unsigned t1 = 0, t2 = mh; t1 < mh; ++t1, ++t2)
            FWHT_2(data[t1], data[t2]);
    }
}

#endif


//------------------------------------------------------------------------------
// Memory Buffer XOR

static void xor_mem(void * LHC_RESTRICT vx, const void * LHC_RESTRICT vy, unsigned bytes)
{
    LHC_M128 * LHC_RESTRICT x16 = reinterpret_cast<LHC_M128 *>(vx);
    const LHC_M128 * LHC_RESTRICT y16 = reinterpret_cast<const LHC_M128 *>(vy);

#if defined(LHC_TARGET_MOBILE)
# if defined(LHC_TRY_NEON)
    // Handle multiples of 64 bytes
    if (CpuHasNeon)
    {
        while (bytes >= 64)
        {
            LHC_M128 x0 = vld1q_u8(x16);
            LHC_M128 x1 = vld1q_u8(x16 + 1);
            LHC_M128 x2 = vld1q_u8(x16 + 2);
            LHC_M128 x3 = vld1q_u8(x16 + 3);
            LHC_M128 y0 = vld1q_u8(y16);
            LHC_M128 y1 = vld1q_u8(y16 + 1);
            LHC_M128 y2 = vld1q_u8(y16 + 2);
            LHC_M128 y3 = vld1q_u8(y16 + 3);

            vst1q_u8(x16,     veorq_u8(x0, y0));
            vst1q_u8(x16 + 1, veorq_u8(x1, y1));
            vst1q_u8(x16 + 2, veorq_u8(x2, y2));
            vst1q_u8(x16 + 3, veorq_u8(x3, y3));

            bytes -= 64, x16 += 4, y16 += 4;
        }

        // Handle multiples of 16 bytes
        while (bytes >= 16)
        {
            LHC_M128 x0 = vld1q_u8(x16);
            LHC_M128 y0 = vld1q_u8(y16);

            vst1q_u8(x16, veorq_u8(x0, y0));

            bytes -= 16, ++x16, ++y16;
        }
    }
    else
# endif // LHC_TRY_NEON
    {
        uint64_t * LHC_RESTRICT x8 = reinterpret_cast<uint64_t *>(x16);
        const uint64_t * LHC_RESTRICT y8 = reinterpret_cast<const uint64_t *>(y16);

        const unsigned count = (unsigned)bytes / 8;
        for (unsigned ii = 0; ii < count; ++ii)
            x8[ii] ^= y8[ii];

        x16 = reinterpret_cast<LHC_M128 *>(x8 + count);
        y16 = reinterpret_cast<const LHC_M128 *>(y8 + count);
    }
#else // LHC_TARGET_MOBILE
# if defined(LHC_TRY_AVX2)
    if (CpuHasAVX2)
    {
        LHC_M256 * LHC_RESTRICT x32 = reinterpret_cast<LHC_M256 *>(x16);
        const LHC_M256 * LHC_RESTRICT y32 = reinterpret_cast<const LHC_M256 *>(y16);

        while (bytes >= 128)
        {
            LHC_M256 x0 = _mm256_loadu_si256(x32);
            LHC_M256 y0 = _mm256_loadu_si256(y32);
            x0 = _mm256_xor_si256(x0, y0);
            LHC_M256 x1 = _mm256_loadu_si256(x32 + 1);
            LHC_M256 y1 = _mm256_loadu_si256(y32 + 1);
            x1 = _mm256_xor_si256(x1, y1);
            LHC_M256 x2 = _mm256_loadu_si256(x32 + 2);
            LHC_M256 y2 = _mm256_loadu_si256(y32 + 2);
            x2 = _mm256_xor_si256(x2, y2);
            LHC_M256 x3 = _mm256_loadu_si256(x32 + 3);
            LHC_M256 y3 = _mm256_loadu_si256(y32 + 3);
            x3 = _mm256_xor_si256(x3, y3);

            _mm256_storeu_si256(x32, x0);
            _mm256_storeu_si256(x32 + 1, x1);
            _mm256_storeu_si256(x32 + 2, x2);
            _mm256_storeu_si256(x32 + 3, x3);

            bytes -= 128, x32 += 4, y32 += 4;
        }

        // Handle multiples of 32 bytes
        while (bytes >= 32)
        {
            // x[i] = x[i] xor y[i]
            _mm256_storeu_si256(x32,
                _mm256_xor_si256(
                    _mm256_loadu_si256(x32),
                    _mm256_loadu_si256(y32)));

            bytes -= 32, ++x32, ++y32;
        }

        x16 = reinterpret_cast<LHC_M128 *>(x32);
        y16 = reinterpret_cast<const LHC_M128 *>(y32);
    }
    else
# endif // LHC_TRY_AVX2
    {
        while (bytes >= 64)
        {
            LHC_M128 x0 = _mm_loadu_si128(x16);
            LHC_M128 y0 = _mm_loadu_si128(y16);
            x0 = _mm_xor_si128(x0, y0);
            LHC_M128 x1 = _mm_loadu_si128(x16 + 1);
            LHC_M128 y1 = _mm_loadu_si128(y16 + 1);
            x1 = _mm_xor_si128(x1, y1);
            LHC_M128 x2 = _mm_loadu_si128(x16 + 2);
            LHC_M128 y2 = _mm_loadu_si128(y16 + 2);
            x2 = _mm_xor_si128(x2, y2);
            LHC_M128 x3 = _mm_loadu_si128(x16 + 3);
            LHC_M128 y3 = _mm_loadu_si128(y16 + 3);
            x3 = _mm_xor_si128(x3, y3);

            _mm_storeu_si128(x16, x0);
            _mm_storeu_si128(x16 + 1, x1);
            _mm_storeu_si128(x16 + 2, x2);
            _mm_storeu_si128(x16 + 3, x3);

            bytes -= 64, x16 += 4, y16 += 4;
        }
    }
#endif // LHC_TARGET_MOBILE

    // Handle multiples of 16 bytes
    while (bytes >= 16)
    {
        // x[i] = x[i] xor y[i]
        _mm_storeu_si128(x16,
            _mm_xor_si128(
                _mm_loadu_si128(x16),
                _mm_loadu_si128(y16)));

        bytes -= 16, ++x16, ++y16;
    }

    uint8_t * LHC_RESTRICT x1 = reinterpret_cast<uint8_t *>(x16);
    const uint8_t * LHC_RESTRICT y1 = reinterpret_cast<const uint8_t *>(y16);

    // Handle a block of 8 bytes
    const unsigned eight = bytes & 8;
    if (eight)
    {
        uint64_t * LHC_RESTRICT x8 = reinterpret_cast<uint64_t *>(x1);
        const uint64_t * LHC_RESTRICT y8 = reinterpret_cast<const uint64_t *>(y1);
        *x8 ^= *y8;
    }

    // Handle a block of 4 bytes
    const unsigned four = bytes & 4;
    if (four)
    {
        uint32_t * LHC_RESTRICT x4 = reinterpret_cast<uint32_t *>(x1 + eight);
        const uint32_t * LHC_RESTRICT y4 = reinterpret_cast<const uint32_t *>(y1 + eight);
        *x4 ^= *y4;
    }

    // Handle final bytes
    const unsigned offset = eight + four;
    switch (bytes & 3)
    {
    case 3: x1[offset + 2] ^= y1[offset + 2];
    case 2: x1[offset + 1] ^= y1[offset + 1];
    case 1: x1[offset] ^= y1[offset];
    default:
        break;
    }
}


//------------------------------------------------------------------------------
// Formal Derivative

// Formal derivative of polynomial in the new basis
static void formal_derivative(GFSymbol* cos, const unsigned size)
{
    for (unsigned i = 1; i < size; ++i)
    {
        const unsigned leng = ((i ^ (i - 1)) + 1) >> 1;

        // If a large number of values are being XORed:
        if (leng >= 8)
            xor_mem(cos + i - leng, cos + i, leng * sizeof(GFSymbol));
        else
            for (unsigned j = i - leng; j < i; j++)
                cos[j] ^= cos[j + leng];
    }

    for (unsigned i = size; i < kFieldSize; i <<= 1)
        xor_mem(cos, cos + i, size * sizeof(GFSymbol));
}


//------------------------------------------------------------------------------
// Fast Fourier Transform

static GFSymbol skewVec[kFieldModulus]; // twisted factors used in FFT

// IFFT in the proposed basis
static void IFLT(GFSymbol* data, const unsigned size, const unsigned index)
{
    for (unsigned depart_no = 1; depart_no < size; depart_no <<= 1)
    {
        for (unsigned j = depart_no; j < size; j += (depart_no << 1))
        {
            // If a large number of values are being XORed:
            if (depart_no >= 8)
                xor_mem(data + j, data + j - depart_no, depart_no * sizeof(GFSymbol));
            else
                for (unsigned i = j - depart_no; i < j; ++i)
                    data[i + depart_no] ^= data[i];

            const GFSymbol skew = skewVec[j + index - 1];

            if (skew != kFieldModulus)
                muladd_mem(data + j - depart_no, data + j, skew, depart_no);
        }
    }
}

// FFT in the proposed basis
static void FLT(GFSymbol* data, const unsigned size, const unsigned index)
{
    for (unsigned depart_no = (size >> 1); depart_no > 0; depart_no >>= 1)
    {
        for (unsigned j = depart_no; j < size; j += (depart_no << 1))
        {
            const GFSymbol skew = skewVec[j + index - 1];

            if (skew != kFieldModulus)
                muladd_mem(data + j - depart_no, data + j, skew, depart_no);

            // If a large number of values are being XORed:
            if (depart_no >= 8)
                xor_mem(data + j, data + j - depart_no, depart_no * sizeof(GFSymbol));
            else
                for (unsigned i = j - depart_no; i < j; ++i)
                    data[i + depart_no] ^= data[i];
        }
    }
}


//------------------------------------------------------------------------------
// FFT Initialization

static GFSymbol B[kFieldSize >> 1];     // factors used in formal derivative
static GFSymbol log_walsh[kFieldSize];  // factors used in the evaluation of the error locator polynomial

// Initialize skewVec[], B[], log_walsh[]
static void InitFieldOperations()
{
	GFSymbol temp[kGFBits - 1];

    for (unsigned i = 1; i < kGFBits; ++i)
        temp[i - 1] = (GFSymbol)((unsigned)1 << i);

    for (unsigned m = 0; m < (kGFBits - 1); ++m)
    {
        const unsigned step = (unsigned)1 << (m + 1);

        skewVec[((unsigned)1 << m) - 1] = 0;

        for (unsigned i = m; i < (kGFBits - 1); ++i)
        {
            const unsigned s = ((unsigned)1 << (i + 1));

            for (unsigned j = ((unsigned)1 << m) - 1; j < s; j += step)
                skewVec[j + s] = skewVec[j] ^ temp[i];
        }

        temp[m] = kFieldModulus - GFLog[mulE(temp[m], GFLog[temp[m] ^ 1])];

        for (unsigned i = m + 1; i < (kGFBits - 1); ++i)
            temp[i] = mulE(temp[i], (GFLog[temp[i] ^ 1] + temp[m]) % kFieldModulus);
    }

    for (unsigned i = 0; i < kFieldSize; ++i)
        skewVec[i] = GFLog[skewVec[i]];

    temp[0] = kFieldModulus - temp[0];

    for (unsigned i = 1; i < (kGFBits - 1); ++i)
        temp[i] = (kFieldModulus - temp[i] + temp[i - 1]) % kFieldModulus;

	B[0] = 0;
    for (unsigned i = 0; i < (kGFBits - 1); ++i)
    {
        const unsigned depart = ((unsigned)1 << i);

        for (unsigned j = 0; j < depart; ++j)
            B[j + depart] = (B[j] + temp[i]) % kFieldModulus;
    }

    for (unsigned i = 0; i < kFieldSize; ++i)
        log_walsh[i] = GFLog[i];

	log_walsh[0] = 0;

    FWHT(log_walsh, kGFBits);
}


//------------------------------------------------------------------------------
// Encoder

// Encoding alg for k/n<0.5: message is a power of two
static void encodeL(GFSymbol* data, const unsigned k, GFSymbol* codeword)
{
	memcpy(codeword, data, sizeof(GFSymbol) * k);

	IFLT(codeword, k, 0);

    for (unsigned i = k; i < kFieldSize; i += k)
    {
        memcpy(&codeword[i], codeword, sizeof(GFSymbol) * k);

        FLT(&codeword[i], k, i);
    }

    memcpy(codeword, data, sizeof(GFSymbol) * k);
}

// Encoding alg for k/n>0.5: parity is a power of two.
// data: message array. parity: parity array. mem: buffer(size>= n-k)
static void encodeH(const GFSymbol* data, const unsigned k, GFSymbol* parity, GFSymbol* mem)
{
	const unsigned t = kFieldSize - k;

    memset(parity, 0, sizeof(GFSymbol) * t);

    for (unsigned i = t; i < kFieldSize; i += t)
    {
        memcpy(mem, &data[i - t], sizeof(GFSymbol) * t);

        IFLT(mem, t, i);

        xor_mem(parity, mem, t * sizeof(GFSymbol));
    }

    FLT(parity, t, 0);
}


//------------------------------------------------------------------------------
// Decoder

static void decode(GFSymbol* codeword, const bool* erasure)
{
    GFSymbol log_walsh2[kFieldSize];

    // Compute the evaluations of the error locator polynomial
    for (unsigned i = 0; i < kFieldSize; ++i)
        log_walsh2[i] = erasure[i] ? 1 : 0;

    FWHT(log_walsh2, kGFBits);

    for (unsigned i = 0; i < kFieldSize; ++i)
        log_walsh2[i] = ((unsigned)log_walsh2[i] * log_walsh[i]) % kFieldModulus;

    FWHT(log_walsh2, kGFBits);

    for (unsigned i = 0; i < kFieldSize; ++i)
        if (erasure[i])
            log_walsh2[i] = kFieldModulus - log_walsh2[i];

    // k2 can be replaced with k
	const unsigned k2 = kFieldSize;

    for (unsigned i = 0; i < kFieldSize; ++i)
        codeword[i] = erasure[i] ? 0 : mulE(codeword[i], log_walsh2[i]);

    IFLT(codeword, kFieldSize, 0);

    // formal derivative
    for (unsigned i = 0; i < kFieldSize; i += 2)
    {
		codeword[i] = mulE(codeword[i], kFieldModulus - B[i >> 1]);
		codeword[i + 1] = mulE(codeword[i + 1], kFieldModulus - B[i >> 1]);
	}

	formal_derivative(codeword, k2);

    for (unsigned i = 0; i < k2; i += 2)
    {
        codeword[i] = mulE(codeword[i], B[i >> 1]);
        codeword[i + 1] = mulE(codeword[i + 1], B[i >> 1]);
    }

	FLT(codeword, k2, 0);

    for (unsigned i = 0; i < k2; ++i)
        codeword[i] = erasure[i] ? mulE(codeword[i], log_walsh2[i]) : 0;
}


//------------------------------------------------------------------------------
// Test Application

void test(unsigned k)
{
	//-----------Generating message----------

    // Message array
	GFSymbol data[kFieldSize] = {0};

    // Filled with random numbers
    for (unsigned i = kFieldSize - k; i < kFieldSize; ++i)
        data[i] = (GFSymbol)rand();


	//---------encoding----------

	GFSymbol codeword[kFieldSize];
	encodeH(&data[kFieldSize - k], k, data, codeword);
	//encodeL(data, k, codeword);

	memcpy(codeword, data, sizeof(GFSymbol) * kFieldSize);


	//--------erasure simulation---------

    // Array indicating erasures
	bool erasure[kFieldSize] = {
        false
    };

    for (unsigned i = k; i < kFieldSize; ++i)
        erasure[i] = true;

    // permuting the erasure array
    for (unsigned i = kFieldSize - 1; i > 0; --i)
    {
        unsigned pos = rand() % (i + 1);

        if (i != pos)
        {
            bool tmp = erasure[i];
            erasure[i] = erasure[pos];
            erasure[pos] = tmp;
        }
    }

    // erasure codeword symbols
    for (unsigned i = 0; i < kFieldSize; ++i)
        if (erasure[i])
            codeword[i] = 0;


	//---------main processing----------
    decode(codeword, erasure);

    // Check the correctness of the result
    for (unsigned i = 0; i < kFieldSize; ++i)
    {
        if (erasure[i] == 1)
        {
            if (data[i] != codeword[i])
            {
                printf("Decoding Error!\n");
                return;
            }
        }
    }

    printf("Decoding is successful!\n");
}


//------------------------------------------------------------------------------
// Entrypoint

int main(int argc, char **argv)
{
    srand((unsigned)time(NULL));

    // Initialize architecture-specific code
    lhc_architecture_init();

    // Fill GFLog table and GFExp table
    InitField();

    // Compute factors used in erasure decoder
    InitFieldOperations();

    for (;;)
    {
        // test(int k), k: message size
        test(kFieldSize / 2);
    }

    return 0;
}
