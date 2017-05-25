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

#include "LeopardCommon.h"

namespace leopard {


//------------------------------------------------------------------------------
// Runtime CPU Architecture Check
//
// Feature checks stolen shamelessly from
// https://github.com/jedisct1/libsodium/blob/master/src/libsodium/sodium/runtime.c

#if defined(HAVE_ANDROID_GETCPUFEATURES)
    #include <cpu-features.h>
#endif

#if defined(LEO_TRY_NEON)
# if defined(IOS) && defined(__ARM_NEON__)
// Requires iPhone 5S or newer
# else
// Remember to add LOCAL_STATIC_LIBRARIES := cpufeatures
bool CpuHasNeon = false; // V6 / V7
bool CpuHasNeon64 = false; // 64-bit
# endif
#endif


#if !defined(LEO_TARGET_MOBILE)

#ifdef _MSC_VER
    #include <intrin.h> // __cpuid
    #pragma warning(disable: 4752) // found Intel(R) Advanced Vector Extensions; consider using /arch:AVX
#endif

#ifdef LEO_TRY_AVX2
bool CpuHasAVX2 = false;
#endif
bool CpuHasSSSE3 = false;

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

#endif // defined(LEO_TARGET_MOBILE)


void InitializeCPUArch()
{
#if defined(LEO_TRY_NEON) && defined(HAVE_ANDROID_GETCPUFEATURES)
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

#if !defined(LEO_TARGET_MOBILE)
    unsigned int cpu_info[4];

    _cpuid(cpu_info, 1);
    CpuHasSSSE3 = ((cpu_info[2] & CPUID_ECX_SSSE3) != 0);

#if defined(LEO_TRY_AVX2)
    _cpuid(cpu_info, 7);
    CpuHasAVX2 = ((cpu_info[1] & CPUID_EBX_AVX2) != 0);
#endif // LEO_TRY_AVX2

#endif // LEO_TARGET_MOBILE
}



// vx[] += vy[] * z
static void muladd_mem(GFSymbol * LEO_RESTRICT vx, const GFSymbol * LEO_RESTRICT vy, GFSymbol z, unsigned symbolCount)
{
    for (unsigned i = 0; i < symbolCount; ++i)
    {
        const GFSymbol a = vy[i];
        if (a == 0)
            continue;

        GFSymbol sum1 = static_cast<GFSymbol>(AddModQ(GFLog[a & 0x0f], z));
        GFSymbol value1 = GFExp[sum1];
        if ((a & 0x0f) == 0)
        {
            value1 = 0;
        }
        GFSymbol sum2 = static_cast<GFSymbol>(AddModQ(GFLog[a & 0xf0], z));
        GFSymbol value2 = GFExp[sum2];
        if ((a & 0xf0) == 0)
        {
            value2 = 0;
        }
        GFSymbol sum3 = static_cast<GFSymbol>(AddModQ(GFLog[a & 0x0f00], z));
        GFSymbol value3 = GFExp[sum3];
        if ((a & 0x0f00) == 0)
        {
            value3 = 0;
        }
        GFSymbol sum4 = static_cast<GFSymbol>(AddModQ(GFLog[a & 0xf000], z));
        GFSymbol value4 = GFExp[sum4];
        if ((a & 0xf000) == 0)
        {
            value4 = 0;
        }

        vx[i] ^= value1;
        vx[i] ^= value2;
        vx[i] ^= value3;
        vx[i] ^= value4;
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
#define LEO_FWHT_OPTIMIZED

typedef GFSymbol fwht_t;

// {a, b} = {a + b, a - b} (Mod Q)
static LEO_FORCE_INLINE void FWHT_2(fwht_t& LEO_RESTRICT a, fwht_t& LEO_RESTRICT b)
{
    const fwht_t sum = AddModQ(a, b);
    const fwht_t dif = SubModQ(a, b);
    a = sum;
    b = dif;
}

/*
    FWHT is a minor slice of the runtime and does not grow with data size,
    but I did attempt a few additional optimizations that failed:

    I've attempted to vectorize (with partial reductions) FWHT_4(data, s),
    which is 70% of the algorithm, but it was slower.  Left in _attic_.

    I've attempted to avoid reductions in all or parts of the FWHT.
    The final modular reduction ends up being slower than the savings.
    Specifically I tried doing it for the whole FWHT and also I tried
    doing it just for the FWHT_2 loop in the main routine, but both
    approaches are slower than partial reductions.

    Replacing word reads with wider reads does speed up the operation, but
    at too high a complexity cost relative to minor perf improvement.
*/

#ifndef LEO_FWHT_OPTIMIZED

// Reference implementation
static void FWHT(fwht_t* data, const unsigned bits)
{
    const unsigned size = (unsigned)(1UL << bits);
    for (unsigned width = 1; width < size; width <<= 1)
        for (unsigned i = 0; i < size; i += (width << 1))
            for (unsigned j = i; j < (width + i); ++j)
                FWHT_2(data[j], data[j + width]);
}

#else

static LEO_FORCE_INLINE void FWHT_4(fwht_t* data)
{
    fwht_t t0 = data[0];
    fwht_t t1 = data[1];
    fwht_t t2 = data[2];
    fwht_t t3 = data[3];
    FWHT_2(t0, t1);
    FWHT_2(t2, t3);
    FWHT_2(t0, t2);
    FWHT_2(t1, t3);
    data[0] = t0;
    data[1] = t1;
    data[2] = t2;
    data[3] = t3;
}

static LEO_FORCE_INLINE void FWHT_4(fwht_t* data, unsigned s)
{
    unsigned x = 0;
    fwht_t t0 = data[x];  x += s;
    fwht_t t1 = data[x];  x += s;
    fwht_t t2 = data[x];  x += s;
    fwht_t t3 = data[x];
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

static inline void FWHT_8(fwht_t* data)
{
    fwht_t t0 = data[0];
    fwht_t t1 = data[1];
    fwht_t t2 = data[2];
    fwht_t t3 = data[3];
    fwht_t t4 = data[4];
    fwht_t t5 = data[5];
    fwht_t t6 = data[6];
    fwht_t t7 = data[7];
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

static inline void FWHT_16(fwht_t* data)
{
    fwht_t t0 = data[0];
    fwht_t t1 = data[1];
    fwht_t t2 = data[2];
    fwht_t t3 = data[3];
    fwht_t t4 = data[4];
    fwht_t t5 = data[5];
    fwht_t t6 = data[6];
    fwht_t t7 = data[7];
    fwht_t t8 = data[8];
    fwht_t t9 = data[9];
    fwht_t t10 = data[10];
    fwht_t t11 = data[11];
    fwht_t t12 = data[12];
    fwht_t t13 = data[13];
    fwht_t t14 = data[14];
    fwht_t t15 = data[15];
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

static void FWHT_SmallData(fwht_t* data, unsigned ldn)
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

// Decimation in time (DIT) version
static void FWHT(fwht_t* data, const unsigned ldn)
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

    for (unsigned ldm = 0; ldm < ldn; ++ldm)
    {
        const unsigned mh = (1UL << ldm);
        for (unsigned t1 = 0, t2 = mh; t1 < mh; ++t1, ++t2)
            FWHT_2(data[t1], data[t2]);
    }
}

#endif


//------------------------------------------------------------------------------
// Memory Buffer XOR

static void xor_mem(void * LEO_RESTRICT vx, const void * LEO_RESTRICT vy, unsigned bytes)
{
    LEO_M128 * LEO_RESTRICT x16 = reinterpret_cast<LEO_M128 *>(vx);
    const LEO_M128 * LEO_RESTRICT y16 = reinterpret_cast<const LEO_M128 *>(vy);

#if defined(LEO_TARGET_MOBILE)
# if defined(LEO_TRY_NEON)
    // Handle multiples of 64 bytes
    if (CpuHasNeon)
    {
        while (bytes >= 64)
        {
            LEO_M128 x0 = vld1q_u8(x16);
            LEO_M128 x1 = vld1q_u8(x16 + 1);
            LEO_M128 x2 = vld1q_u8(x16 + 2);
            LEO_M128 x3 = vld1q_u8(x16 + 3);
            LEO_M128 y0 = vld1q_u8(y16);
            LEO_M128 y1 = vld1q_u8(y16 + 1);
            LEO_M128 y2 = vld1q_u8(y16 + 2);
            LEO_M128 y3 = vld1q_u8(y16 + 3);

            vst1q_u8(x16,     veorq_u8(x0, y0));
            vst1q_u8(x16 + 1, veorq_u8(x1, y1));
            vst1q_u8(x16 + 2, veorq_u8(x2, y2));
            vst1q_u8(x16 + 3, veorq_u8(x3, y3));

            bytes -= 64, x16 += 4, y16 += 4;
        }

        // Handle multiples of 16 bytes
        while (bytes >= 16)
        {
            LEO_M128 x0 = vld1q_u8(x16);
            LEO_M128 y0 = vld1q_u8(y16);

            vst1q_u8(x16, veorq_u8(x0, y0));

            bytes -= 16, ++x16, ++y16;
        }
    }
    else
# endif // LEO_TRY_NEON
    {
        uint64_t * LEO_RESTRICT x8 = reinterpret_cast<uint64_t *>(x16);
        const uint64_t * LEO_RESTRICT y8 = reinterpret_cast<const uint64_t *>(y16);

        const unsigned count = (unsigned)bytes / 8;
        for (unsigned ii = 0; ii < count; ++ii)
            x8[ii] ^= y8[ii];

        x16 = reinterpret_cast<LEO_M128 *>(x8 + count);
        y16 = reinterpret_cast<const LEO_M128 *>(y8 + count);
    }
#else // LEO_TARGET_MOBILE
# if defined(LEO_TRY_AVX2)
    if (CpuHasAVX2)
    {
        LEO_M256 * LEO_RESTRICT x32 = reinterpret_cast<LEO_M256 *>(x16);
        const LEO_M256 * LEO_RESTRICT y32 = reinterpret_cast<const LEO_M256 *>(y16);

        while (bytes >= 128)
        {
            LEO_M256 x0 = _mm256_loadu_si256(x32);
            LEO_M256 y0 = _mm256_loadu_si256(y32);
            x0 = _mm256_xor_si256(x0, y0);
            LEO_M256 x1 = _mm256_loadu_si256(x32 + 1);
            LEO_M256 y1 = _mm256_loadu_si256(y32 + 1);
            x1 = _mm256_xor_si256(x1, y1);
            LEO_M256 x2 = _mm256_loadu_si256(x32 + 2);
            LEO_M256 y2 = _mm256_loadu_si256(y32 + 2);
            x2 = _mm256_xor_si256(x2, y2);
            LEO_M256 x3 = _mm256_loadu_si256(x32 + 3);
            LEO_M256 y3 = _mm256_loadu_si256(y32 + 3);
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

        x16 = reinterpret_cast<LEO_M128 *>(x32);
        y16 = reinterpret_cast<const LEO_M128 *>(y32);
    }
    else
# endif // LEO_TRY_AVX2
    {
        while (bytes >= 64)
        {
            LEO_M128 x0 = _mm_loadu_si128(x16);
            LEO_M128 y0 = _mm_loadu_si128(y16);
            x0 = _mm_xor_si128(x0, y0);
            LEO_M128 x1 = _mm_loadu_si128(x16 + 1);
            LEO_M128 y1 = _mm_loadu_si128(y16 + 1);
            x1 = _mm_xor_si128(x1, y1);
            LEO_M128 x2 = _mm_loadu_si128(x16 + 2);
            LEO_M128 y2 = _mm_loadu_si128(y16 + 2);
            x2 = _mm_xor_si128(x2, y2);
            LEO_M128 x3 = _mm_loadu_si128(x16 + 3);
            LEO_M128 y3 = _mm_loadu_si128(y16 + 3);
            x3 = _mm_xor_si128(x3, y3);

            _mm_storeu_si128(x16, x0);
            _mm_storeu_si128(x16 + 1, x1);
            _mm_storeu_si128(x16 + 2, x2);
            _mm_storeu_si128(x16 + 3, x3);

            bytes -= 64, x16 += 4, y16 += 4;
        }
    }
#endif // LEO_TARGET_MOBILE

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

    uint8_t * LEO_RESTRICT x1 = reinterpret_cast<uint8_t *>(x16);
    const uint8_t * LEO_RESTRICT y1 = reinterpret_cast<const uint8_t *>(y16);

    // Handle a block of 8 bytes
    const unsigned eight = bytes & 8;
    if (eight)
    {
        uint64_t * LEO_RESTRICT x8 = reinterpret_cast<uint64_t *>(x1);
        const uint64_t * LEO_RESTRICT y8 = reinterpret_cast<const uint64_t *>(y1);
        *x8 ^= *y8;
    }

    // Handle a block of 4 bytes
    const unsigned four = bytes & 4;
    if (four)
    {
        uint32_t * LEO_RESTRICT x4 = reinterpret_cast<uint32_t *>(x1 + eight);
        const uint32_t * LEO_RESTRICT y4 = reinterpret_cast<const uint32_t *>(y1 + eight);
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
static fwht_t log_walsh[kFieldSize];  // factors used in the evaluation of the error locator polynomial

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

static void decode(GFSymbol* codeword, unsigned k, const bool* erasure)
{
    fwht_t log_walsh2[kFieldSize];

    // Compute the evaluations of the error locator polynomial
    for (unsigned i = 0; i < kFieldSize; ++i)
        log_walsh2[i] = erasure[i] ? 1 : 0;

    FWHT(log_walsh2, kGFBits);

    for (unsigned i = 0; i < kFieldSize; ++i)
        log_walsh2[i] = ((unsigned)log_walsh2[i] * (unsigned)log_walsh[i]) % kFieldModulus;

    FWHT(log_walsh2, kGFBits);

    // k2 can be replaced with k
    const unsigned k2 = kFieldSize;
    //const unsigned k2 = k; // cannot actually be replaced with k.  what else need to change?

    for (unsigned i = 0; i < kFieldSize; ++i)
    {
        if (erasure[i])
        {
            codeword[i] = 0;
        }
        else
        {
            codeword[i] = mulE(codeword[i], log_walsh2[i]);
        }
    }

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
    {
        if (erasure[i])
        {
            codeword[i] = mulE(codeword[i], kFieldModulus - log_walsh2[i]);
        }
    }
}


//------------------------------------------------------------------------------
// Test Application

void test(unsigned k, unsigned seed)
{
    srand(seed);

    //-----------Generating message----------

    // Message array
    GFSymbol data[kFieldSize] = {0};

    // Filled with random numbers
    for (unsigned i = kFieldSize - k; i < kFieldSize; ++i)
        data[i] = (GFSymbol)rand();


    //---------encoding----------

    GFSymbol codeword[kFieldSize];
    encodeH(&data[kFieldSize - k], k, data, codeword);
    //encodeL(data, k, codeword); // does not seem to work with any input?  what else needs to change?

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
    decode(codeword, k, erasure);

    // Check the correctness of the result
    for (unsigned i = 0; i < kFieldSize; ++i)
    {
        if (erasure[i] == 1)
        {
            if (data[i] != codeword[i])
            {
                printf("Decoding Error with seed = %d!\n", seed);
                LEO_DEBUG_BREAK;
                return;
            }
        }
    }

    //printf("Decoding is successful!\n");
}


//------------------------------------------------------------------------------
// Entrypoint

int main(int argc, char **argv)
{
    // Initialize architecture-specific code
    leo_architecture_init();

    // Fill GFLog table and GFExp table
    InitField();

    // Compute factors used in erasure decoder
    InitFieldOperations();

    unsigned seed = (unsigned)time(NULL);
    for (;;)
    {
        // test(int k), k: message size
        /*
            EncodeH works for kFieldSize / 2 and kFieldSize * 3 / 4, etc,
            s.t. the number of recovery pieces is a power of two
        */
        test(kFieldSize / 2, seed);

        ++seed;
    }

    return 0;
}


} // namespace leopard
