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
#include "Leopard2Backend.h"

#include <thread>

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
#define CPUID_ECX_OSXSAVE 0x08000000
#define CPUID_ECX_AVX     0x10000000

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

#if defined(LEO_TRY_AVX2)
static uint64_t ReadXCR0()
{
#if defined(_MSC_VER) && (defined(_M_X64) || defined(_M_AMD64) || defined(_M_IX86))
    return static_cast<uint64_t>(_xgetbv(0));
#elif defined(__i386__) || defined(__x86_64__)
    unsigned int eax, edx;
    __asm__ __volatile__ ("xgetbv" : "=a" (eax), "=d" (edx) : "c" (0));
    return (static_cast<uint64_t>(edx) << 32) | eax;
#else
    return 0;
#endif
}
#endif

#elif defined(LEO_USE_SSE2NEON)
bool CpuHasSSSE3 = true;
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
    unsigned int max_basic_leaf = 0;

    _cpuid(cpu_info, 0);
    max_basic_leaf = cpu_info[0];

#if defined(LEO_TRY_AVX2)
    bool ymm_state_enabled = false;
#endif
    if (max_basic_leaf >= 1)
    {
        _cpuid(cpu_info, 1);
        CpuHasSSSE3 = ((cpu_info[2] & CPUID_ECX_SSSE3) != 0);

#if defined(LEO_TRY_AVX2)
        const unsigned int avx_os_mask = CPUID_ECX_AVX | CPUID_ECX_OSXSAVE;
        if ((cpu_info[2] & avx_os_mask) == avx_os_mask)
        {
            // XMM state (bit 1) and YMM state (bit 2) must both be managed by
            // the OS before any AVX/AVX2 instruction can execute safely.
            ymm_state_enabled = (ReadXCR0() & 0x6) == 0x6;
        }
#endif
    }

#if defined(LEO_TRY_AVX2)
    if (ymm_state_enabled && max_basic_leaf >= 7)
    {
        _cpuid(cpu_info, 7);
        CpuHasAVX2 = ((cpu_info[1] & CPUID_EBX_AVX2) != 0);
    }
#endif // LEO_TRY_AVX2

#if !defined(LEO_USE_SSSE3_OPT) || !defined(LEO_TRY_SSSE3)
    CpuHasSSSE3 = false;
#endif
#if defined(LEO_TRY_AVX2) && !defined(LEO_USE_AVX2_OPT)
    CpuHasAVX2 = false;
#endif

#endif // LEO_TARGET_MOBILE
}


namespace backend {

leo2_backend ExecutionBackend()
{
    if (SelectedBackend() == LEO2_BACKEND_AUTO)
        return LEO2_BACKEND_AUTO;
#if defined(LEO_TRY_NEON)
    if (CpuHasNeon)
        return LEO2_BACKEND_NEON;
#endif
#if defined(LEO_USE_SSE2NEON)
    if (CpuHasSSSE3)
        return LEO2_BACKEND_NEON;
#endif
    return SelectedBackend();
}

} // namespace backend


//------------------------------------------------------------------------------
// XOR Memory

void xor_mem_baseline(
    void * LEO_RESTRICT vx, const void * LEO_RESTRICT vy,
    uint64_t bytes)
{
#if defined(LEO_TRY_AVX2)
    if (CpuHasAVX2)
    {
        LEO_M256 * LEO_RESTRICT x32 = reinterpret_cast<LEO_M256 *>(vx);
        const LEO_M256 * LEO_RESTRICT y32 = reinterpret_cast<const LEO_M256 *>(vy);
        while (bytes >= 128)
        {
            const LEO_M256 x0 = _mm256_xor_si256(_mm256_loadu_si256(x32),     _mm256_loadu_si256(y32));
            const LEO_M256 x1 = _mm256_xor_si256(_mm256_loadu_si256(x32 + 1), _mm256_loadu_si256(y32 + 1));
            const LEO_M256 x2 = _mm256_xor_si256(_mm256_loadu_si256(x32 + 2), _mm256_loadu_si256(y32 + 2));
            const LEO_M256 x3 = _mm256_xor_si256(_mm256_loadu_si256(x32 + 3), _mm256_loadu_si256(y32 + 3));
            _mm256_storeu_si256(x32, x0);
            _mm256_storeu_si256(x32 + 1, x1);
            _mm256_storeu_si256(x32 + 2, x2);
            _mm256_storeu_si256(x32 + 3, x3);
            x32 += 4, y32 += 4;
            bytes -= 128;
        };
        if (bytes > 0)
        {
            const LEO_M256 x0 = _mm256_xor_si256(_mm256_loadu_si256(x32),     _mm256_loadu_si256(y32));
            const LEO_M256 x1 = _mm256_xor_si256(_mm256_loadu_si256(x32 + 1), _mm256_loadu_si256(y32 + 1));
            _mm256_storeu_si256(x32, x0);
            _mm256_storeu_si256(x32 + 1, x1);
        }
        return;
    }
#endif // LEO_TRY_AVX2

    uint8_t* LEO_RESTRICT x_bytes = reinterpret_cast<uint8_t*>(vx);
    const uint8_t* LEO_RESTRICT y_bytes = reinterpret_cast<const uint8_t*>(vy);
    while (bytes >= 64)
    {
        const LEO_M128 x0 = _mm_xor_si128(
            _mm_loadu_si128(reinterpret_cast<const LEO_M128*>(x_bytes)),
            _mm_loadu_si128(reinterpret_cast<const LEO_M128*>(y_bytes)));
        const LEO_M128 x1 = _mm_xor_si128(
            _mm_loadu_si128(reinterpret_cast<const LEO_M128*>(x_bytes + 16)),
            _mm_loadu_si128(reinterpret_cast<const LEO_M128*>(y_bytes + 16)));
        const LEO_M128 x2 = _mm_xor_si128(
            _mm_loadu_si128(reinterpret_cast<const LEO_M128*>(x_bytes + 32)),
            _mm_loadu_si128(reinterpret_cast<const LEO_M128*>(y_bytes + 32)));
        const LEO_M128 x3 = _mm_xor_si128(
            _mm_loadu_si128(reinterpret_cast<const LEO_M128*>(x_bytes + 48)),
            _mm_loadu_si128(reinterpret_cast<const LEO_M128*>(y_bytes + 48)));
        _mm_storeu_si128(reinterpret_cast<LEO_M128*>(x_bytes), x0);
        _mm_storeu_si128(reinterpret_cast<LEO_M128*>(x_bytes + 16), x1);
        _mm_storeu_si128(reinterpret_cast<LEO_M128*>(x_bytes + 32), x2);
        _mm_storeu_si128(reinterpret_cast<LEO_M128*>(x_bytes + 48), x3);
        x_bytes += 64;
        y_bytes += 64;
        bytes -= 64;
    }
    while (bytes >= 16)
    {
        const LEO_M128 result = _mm_xor_si128(
            _mm_loadu_si128(reinterpret_cast<const LEO_M128*>(x_bytes)),
            _mm_loadu_si128(reinterpret_cast<const LEO_M128*>(y_bytes)));
        _mm_storeu_si128(reinterpret_cast<LEO_M128*>(x_bytes), result);
        x_bytes += 16;
        y_bytes += 16;
        bytes -= 16;
    }
    while (bytes-- != 0)
        *x_bytes++ ^= *y_bytes++;
}

void xor_mem(
    void * LEO_RESTRICT vx, const void * LEO_RESTRICT vy,
    uint64_t bytes)
{
    if (backend::SelectedBackend() != LEO2_BACKEND_AUTO)
        backend::GetDefaultOps().xor_memory(vx, vy, bytes);
    else
        xor_mem_baseline(vx, vy, bytes);
}

void xor_mem(
    const backend::Ops& ops,
    void * LEO_RESTRICT vx, const void * LEO_RESTRICT vy,
    uint64_t bytes)
{
    ops.xor_memory(vx, vy, bytes);
}

#ifdef LEO_M1_OPT

void xor_mem_2to1(
    void * LEO_RESTRICT x,
    const void * LEO_RESTRICT y,
    const void * LEO_RESTRICT z,
    uint64_t bytes)
{
#if defined(LEO_TRY_AVX2)
    if (CpuHasAVX2)
    {
        LEO_M256 * LEO_RESTRICT x32 = reinterpret_cast<LEO_M256 *>(x);
        const LEO_M256 * LEO_RESTRICT y32 = reinterpret_cast<const LEO_M256 *>(y);
        const LEO_M256 * LEO_RESTRICT z32 = reinterpret_cast<const LEO_M256 *>(z);
        while (bytes >= 128)
        {
            LEO_M256 x0 = _mm256_xor_si256(_mm256_loadu_si256(x32), _mm256_loadu_si256(y32));
            x0 = _mm256_xor_si256(x0, _mm256_loadu_si256(z32));
            LEO_M256 x1 = _mm256_xor_si256(_mm256_loadu_si256(x32 + 1), _mm256_loadu_si256(y32 + 1));
            x1 = _mm256_xor_si256(x1, _mm256_loadu_si256(z32 + 1));
            LEO_M256 x2 = _mm256_xor_si256(_mm256_loadu_si256(x32 + 2), _mm256_loadu_si256(y32 + 2));
            x2 = _mm256_xor_si256(x2, _mm256_loadu_si256(z32 + 2));
            LEO_M256 x3 = _mm256_xor_si256(_mm256_loadu_si256(x32 + 3), _mm256_loadu_si256(y32 + 3));
            x3 = _mm256_xor_si256(x3, _mm256_loadu_si256(z32 + 3));
            _mm256_storeu_si256(x32, x0);
            _mm256_storeu_si256(x32 + 1, x1);
            _mm256_storeu_si256(x32 + 2, x2);
            _mm256_storeu_si256(x32 + 3, x3);
            x32 += 4, y32 += 4, z32 += 4;
            bytes -= 128;
        };

        if (bytes > 0)
        {
            LEO_M256 x0 = _mm256_xor_si256(_mm256_loadu_si256(x32),     _mm256_loadu_si256(y32));
            x0 = _mm256_xor_si256(x0, _mm256_loadu_si256(z32));
            LEO_M256 x1 = _mm256_xor_si256(_mm256_loadu_si256(x32 + 1), _mm256_loadu_si256(y32 + 1));
            x1 = _mm256_xor_si256(x1, _mm256_loadu_si256(z32 + 1));
            _mm256_storeu_si256(x32, x0);
            _mm256_storeu_si256(x32 + 1, x1);
        }

        return;
    }
#endif // LEO_TRY_AVX2

    LEO_M128 * LEO_RESTRICT x16 = reinterpret_cast<LEO_M128 *>(x);
    const LEO_M128 * LEO_RESTRICT y16 = reinterpret_cast<const LEO_M128 *>(y);
    const LEO_M128 * LEO_RESTRICT z16 = reinterpret_cast<const LEO_M128 *>(z);
    do
    {
        LEO_M128 x0 = _mm_xor_si128(_mm_loadu_si128(x16), _mm_loadu_si128(y16));
        x0 = _mm_xor_si128(x0, _mm_loadu_si128(z16));
        LEO_M128 x1 = _mm_xor_si128(_mm_loadu_si128(x16 + 1), _mm_loadu_si128(y16 + 1));
        x1 = _mm_xor_si128(x1, _mm_loadu_si128(z16 + 1));
        LEO_M128 x2 = _mm_xor_si128(_mm_loadu_si128(x16 + 2), _mm_loadu_si128(y16 + 2));
        x2 = _mm_xor_si128(x2, _mm_loadu_si128(z16 + 2));
        LEO_M128 x3 = _mm_xor_si128(_mm_loadu_si128(x16 + 3), _mm_loadu_si128(y16 + 3));
        x3 = _mm_xor_si128(x3, _mm_loadu_si128(z16 + 3));
        _mm_storeu_si128(x16, x0);
        _mm_storeu_si128(x16 + 1, x1);
        _mm_storeu_si128(x16 + 2, x2);
        _mm_storeu_si128(x16 + 3, x3);
        x16 += 4, y16 += 4, z16 += 4;
        bytes -= 64;
    } while (bytes > 0);
}

void xor_mem_2to1(
    const backend::Ops& ops,
    void * LEO_RESTRICT x,
    const void * LEO_RESTRICT y,
    const void * LEO_RESTRICT z,
    uint64_t bytes)
{
    // Preserve the fused one-pass legacy path for AUTO/process-default calls.
    // An explicitly lower (or test-tracing) table must remain authoritative,
    // even though that requires two XOR passes until Ops grows a fused slot.
    if (&ops == &backend::GetDefaultOps())
        xor_mem_2to1(x, y, z, bytes);
    else
    {
        ops.xor_memory(x, y, bytes);
        ops.xor_memory(x, z, bytes);
    }
}

#endif // LEO_M1_OPT

#ifdef LEO_USE_VECTOR4_OPT

void xor_mem4(
    void * LEO_RESTRICT vx_0, const void * LEO_RESTRICT vy_0,
    void * LEO_RESTRICT vx_1, const void * LEO_RESTRICT vy_1,
    void * LEO_RESTRICT vx_2, const void * LEO_RESTRICT vy_2,
    void * LEO_RESTRICT vx_3, const void * LEO_RESTRICT vy_3,
    uint64_t bytes)
{
#if defined(LEO_TRY_AVX2)
    if (CpuHasAVX2)
    {
        LEO_M256 * LEO_RESTRICT       x32_0 = reinterpret_cast<LEO_M256 *>      (vx_0);
        const LEO_M256 * LEO_RESTRICT y32_0 = reinterpret_cast<const LEO_M256 *>(vy_0);
        LEO_M256 * LEO_RESTRICT       x32_1 = reinterpret_cast<LEO_M256 *>      (vx_1);
        const LEO_M256 * LEO_RESTRICT y32_1 = reinterpret_cast<const LEO_M256 *>(vy_1);
        LEO_M256 * LEO_RESTRICT       x32_2 = reinterpret_cast<LEO_M256 *>      (vx_2);
        const LEO_M256 * LEO_RESTRICT y32_2 = reinterpret_cast<const LEO_M256 *>(vy_2);
        LEO_M256 * LEO_RESTRICT       x32_3 = reinterpret_cast<LEO_M256 *>      (vx_3);
        const LEO_M256 * LEO_RESTRICT y32_3 = reinterpret_cast<const LEO_M256 *>(vy_3);
        while (bytes >= 128)
        {
            const LEO_M256 x0_0 = _mm256_xor_si256(_mm256_loadu_si256(x32_0),     _mm256_loadu_si256(y32_0));
            const LEO_M256 x1_0 = _mm256_xor_si256(_mm256_loadu_si256(x32_0 + 1), _mm256_loadu_si256(y32_0 + 1));
            const LEO_M256 x2_0 = _mm256_xor_si256(_mm256_loadu_si256(x32_0 + 2), _mm256_loadu_si256(y32_0 + 2));
            const LEO_M256 x3_0 = _mm256_xor_si256(_mm256_loadu_si256(x32_0 + 3), _mm256_loadu_si256(y32_0 + 3));
            _mm256_storeu_si256(x32_0, x0_0);
            _mm256_storeu_si256(x32_0 + 1, x1_0);
            _mm256_storeu_si256(x32_0 + 2, x2_0);
            _mm256_storeu_si256(x32_0 + 3, x3_0);
            x32_0 += 4, y32_0 += 4;
            const LEO_M256 x0_1 = _mm256_xor_si256(_mm256_loadu_si256(x32_1),     _mm256_loadu_si256(y32_1));
            const LEO_M256 x1_1 = _mm256_xor_si256(_mm256_loadu_si256(x32_1 + 1), _mm256_loadu_si256(y32_1 + 1));
            const LEO_M256 x2_1 = _mm256_xor_si256(_mm256_loadu_si256(x32_1 + 2), _mm256_loadu_si256(y32_1 + 2));
            const LEO_M256 x3_1 = _mm256_xor_si256(_mm256_loadu_si256(x32_1 + 3), _mm256_loadu_si256(y32_1 + 3));
            _mm256_storeu_si256(x32_1, x0_1);
            _mm256_storeu_si256(x32_1 + 1, x1_1);
            _mm256_storeu_si256(x32_1 + 2, x2_1);
            _mm256_storeu_si256(x32_1 + 3, x3_1);
            x32_1 += 4, y32_1 += 4;
            const LEO_M256 x0_2 = _mm256_xor_si256(_mm256_loadu_si256(x32_2),     _mm256_loadu_si256(y32_2));
            const LEO_M256 x1_2 = _mm256_xor_si256(_mm256_loadu_si256(x32_2 + 1), _mm256_loadu_si256(y32_2 + 1));
            const LEO_M256 x2_2 = _mm256_xor_si256(_mm256_loadu_si256(x32_2 + 2), _mm256_loadu_si256(y32_2 + 2));
            const LEO_M256 x3_2 = _mm256_xor_si256(_mm256_loadu_si256(x32_2 + 3), _mm256_loadu_si256(y32_2 + 3));
            _mm256_storeu_si256(x32_2, x0_2);
            _mm256_storeu_si256(x32_2 + 1, x1_2);
            _mm256_storeu_si256(x32_2 + 2, x2_2);
            _mm256_storeu_si256(x32_2 + 3, x3_2);
            x32_2 += 4, y32_2 += 4;
            const LEO_M256 x0_3 = _mm256_xor_si256(_mm256_loadu_si256(x32_3),     _mm256_loadu_si256(y32_3));
            const LEO_M256 x1_3 = _mm256_xor_si256(_mm256_loadu_si256(x32_3 + 1), _mm256_loadu_si256(y32_3 + 1));
            const LEO_M256 x2_3 = _mm256_xor_si256(_mm256_loadu_si256(x32_3 + 2), _mm256_loadu_si256(y32_3 + 2));
            const LEO_M256 x3_3 = _mm256_xor_si256(_mm256_loadu_si256(x32_3 + 3), _mm256_loadu_si256(y32_3 + 3));
            _mm256_storeu_si256(x32_3,     x0_3);
            _mm256_storeu_si256(x32_3 + 1, x1_3);
            _mm256_storeu_si256(x32_3 + 2, x2_3);
            _mm256_storeu_si256(x32_3 + 3, x3_3);
            x32_3 += 4, y32_3 += 4;
            bytes -= 128;
        }
        if (bytes > 0)
        {
            const LEO_M256 x0_0 = _mm256_xor_si256(_mm256_loadu_si256(x32_0),     _mm256_loadu_si256(y32_0));
            const LEO_M256 x1_0 = _mm256_xor_si256(_mm256_loadu_si256(x32_0 + 1), _mm256_loadu_si256(y32_0 + 1));
            const LEO_M256 x0_1 = _mm256_xor_si256(_mm256_loadu_si256(x32_1),     _mm256_loadu_si256(y32_1));
            const LEO_M256 x1_1 = _mm256_xor_si256(_mm256_loadu_si256(x32_1 + 1), _mm256_loadu_si256(y32_1 + 1));
            _mm256_storeu_si256(x32_0, x0_0);
            _mm256_storeu_si256(x32_0 + 1, x1_0);
            _mm256_storeu_si256(x32_1, x0_1);
            _mm256_storeu_si256(x32_1 + 1, x1_1);
            const LEO_M256 x0_2 = _mm256_xor_si256(_mm256_loadu_si256(x32_2),     _mm256_loadu_si256(y32_2));
            const LEO_M256 x1_2 = _mm256_xor_si256(_mm256_loadu_si256(x32_2 + 1), _mm256_loadu_si256(y32_2 + 1));
            const LEO_M256 x0_3 = _mm256_xor_si256(_mm256_loadu_si256(x32_3),     _mm256_loadu_si256(y32_3));
            const LEO_M256 x1_3 = _mm256_xor_si256(_mm256_loadu_si256(x32_3 + 1), _mm256_loadu_si256(y32_3 + 1));
            _mm256_storeu_si256(x32_2,     x0_2);
            _mm256_storeu_si256(x32_2 + 1, x1_2);
            _mm256_storeu_si256(x32_3,     x0_3);
            _mm256_storeu_si256(x32_3 + 1, x1_3);
        }
        return;
    }
#endif // LEO_TRY_AVX2
    LEO_M128 * LEO_RESTRICT       x16_0 = reinterpret_cast<LEO_M128 *>      (vx_0);
    const LEO_M128 * LEO_RESTRICT y16_0 = reinterpret_cast<const LEO_M128 *>(vy_0);
    LEO_M128 * LEO_RESTRICT       x16_1 = reinterpret_cast<LEO_M128 *>      (vx_1);
    const LEO_M128 * LEO_RESTRICT y16_1 = reinterpret_cast<const LEO_M128 *>(vy_1);
    LEO_M128 * LEO_RESTRICT       x16_2 = reinterpret_cast<LEO_M128 *>      (vx_2);
    const LEO_M128 * LEO_RESTRICT y16_2 = reinterpret_cast<const LEO_M128 *>(vy_2);
    LEO_M128 * LEO_RESTRICT       x16_3 = reinterpret_cast<LEO_M128 *>      (vx_3);
    const LEO_M128 * LEO_RESTRICT y16_3 = reinterpret_cast<const LEO_M128 *>(vy_3);
    do
    {
        const LEO_M128 x0_0 = _mm_xor_si128(_mm_loadu_si128(x16_0),     _mm_loadu_si128(y16_0));
        const LEO_M128 x1_0 = _mm_xor_si128(_mm_loadu_si128(x16_0 + 1), _mm_loadu_si128(y16_0 + 1));
        const LEO_M128 x2_0 = _mm_xor_si128(_mm_loadu_si128(x16_0 + 2), _mm_loadu_si128(y16_0 + 2));
        const LEO_M128 x3_0 = _mm_xor_si128(_mm_loadu_si128(x16_0 + 3), _mm_loadu_si128(y16_0 + 3));
        _mm_storeu_si128(x16_0, x0_0);
        _mm_storeu_si128(x16_0 + 1, x1_0);
        _mm_storeu_si128(x16_0 + 2, x2_0);
        _mm_storeu_si128(x16_0 + 3, x3_0);
        x16_0 += 4, y16_0 += 4;
        const LEO_M128 x0_1 = _mm_xor_si128(_mm_loadu_si128(x16_1),     _mm_loadu_si128(y16_1));
        const LEO_M128 x1_1 = _mm_xor_si128(_mm_loadu_si128(x16_1 + 1), _mm_loadu_si128(y16_1 + 1));
        const LEO_M128 x2_1 = _mm_xor_si128(_mm_loadu_si128(x16_1 + 2), _mm_loadu_si128(y16_1 + 2));
        const LEO_M128 x3_1 = _mm_xor_si128(_mm_loadu_si128(x16_1 + 3), _mm_loadu_si128(y16_1 + 3));
        _mm_storeu_si128(x16_1, x0_1);
        _mm_storeu_si128(x16_1 + 1, x1_1);
        _mm_storeu_si128(x16_1 + 2, x2_1);
        _mm_storeu_si128(x16_1 + 3, x3_1);
        x16_1 += 4, y16_1 += 4;
        const LEO_M128 x0_2 = _mm_xor_si128(_mm_loadu_si128(x16_2),     _mm_loadu_si128(y16_2));
        const LEO_M128 x1_2 = _mm_xor_si128(_mm_loadu_si128(x16_2 + 1), _mm_loadu_si128(y16_2 + 1));
        const LEO_M128 x2_2 = _mm_xor_si128(_mm_loadu_si128(x16_2 + 2), _mm_loadu_si128(y16_2 + 2));
        const LEO_M128 x3_2 = _mm_xor_si128(_mm_loadu_si128(x16_2 + 3), _mm_loadu_si128(y16_2 + 3));
        _mm_storeu_si128(x16_2, x0_2);
        _mm_storeu_si128(x16_2 + 1, x1_2);
        _mm_storeu_si128(x16_2 + 2, x2_2);
        _mm_storeu_si128(x16_2 + 3, x3_2);
        x16_2 += 4, y16_2 += 4;
        const LEO_M128 x0_3 = _mm_xor_si128(_mm_loadu_si128(x16_3),     _mm_loadu_si128(y16_3));
        const LEO_M128 x1_3 = _mm_xor_si128(_mm_loadu_si128(x16_3 + 1), _mm_loadu_si128(y16_3 + 1));
        const LEO_M128 x2_3 = _mm_xor_si128(_mm_loadu_si128(x16_3 + 2), _mm_loadu_si128(y16_3 + 2));
        const LEO_M128 x3_3 = _mm_xor_si128(_mm_loadu_si128(x16_3 + 3), _mm_loadu_si128(y16_3 + 3));
        _mm_storeu_si128(x16_3,     x0_3);
        _mm_storeu_si128(x16_3 + 1, x1_3);
        _mm_storeu_si128(x16_3 + 2, x2_3);
        _mm_storeu_si128(x16_3 + 3, x3_3);
        x16_3 += 4, y16_3 += 4;
        bytes -= 64;
    } while (bytes > 0);
}

#endif // LEO_USE_VECTOR4_OPT

void VectorXOR_Threads(
    const uint64_t bytes,
    unsigned count,
    void** x,
    void** y)
{
    VectorXOR_Threads(backend::GetDefaultOps(), bytes, count, x, y);
}

void VectorXOR_Threads(
    const backend::Ops& ops,
    const uint64_t bytes,
    unsigned count,
    void** x,
    void** y)
{
#ifdef LEO_USE_VECTOR4_OPT
    if (count >= 4)
    {
        int i_end = count - 4;
#pragma omp parallel for
        for (int i = 0; i <= i_end; i += 4)
        {
            xor_mem4(
                x[i + 0], y[i + 0],
                x[i + 1], y[i + 1],
                x[i + 2], y[i + 2],
                x[i + 3], y[i + 3],
                bytes);
        }
        count %= 4;
        i_end -= count;
        x += i_end;
        y += i_end;
    }
#endif // LEO_USE_VECTOR4_OPT

    for (unsigned i = 0; i < count; ++i)
        xor_mem(ops, x[i], y[i], bytes);
}

void VectorXOR(
    const uint64_t bytes,
    unsigned count,
    void** x,
    void** y)
{
    VectorXOR(backend::GetDefaultOps(), bytes, count, x, y);
}

void VectorXOR(
    const backend::Ops& ops,
    const uint64_t bytes,
    unsigned count,
    void** x,
    void** y)
{
#ifdef LEO_USE_VECTOR4_OPT
    if (count >= 4)
    {
        int i_end = count - 4;
        for (int i = 0; i <= i_end; i += 4)
        {
            xor_mem4(
                x[i + 0], y[i + 0],
                x[i + 1], y[i + 1],
                x[i + 2], y[i + 2],
                x[i + 3], y[i + 3],
                bytes);
        }
        count %= 4;
        i_end -= count;
        x += i_end;
        y += i_end;
    }
#endif // LEO_USE_VECTOR4_OPT

    for (unsigned i = 0; i < count; ++i)
        xor_mem(ops, x[i], y[i], bytes);
}


} // namespace leopard
