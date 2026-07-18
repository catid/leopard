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

#include "benchmark_ff8xor_xor3.h"

#include <errno.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#if defined(__linux__)
# include <sched.h>
# include <unistd.h>
#endif

#if defined(_MSC_VER)
# include <intrin.h>
#endif


namespace {


struct CpuFeatures
{
    bool avx512f;
    bool avx512vl;
    bool avx2;
    bool rdpru;
};

static void Cpuid(
    uint32_t leaf,
    uint32_t subleaf,
    uint32_t& eax,
    uint32_t& ebx,
    uint32_t& ecx,
    uint32_t& edx)
{
#if defined(_MSC_VER) && (defined(_M_X64) || defined(_M_IX86))
    int registers[4];
    __cpuidex(registers, static_cast<int>(leaf), static_cast<int>(subleaf));
    eax = static_cast<uint32_t>(registers[0]);
    ebx = static_cast<uint32_t>(registers[1]);
    ecx = static_cast<uint32_t>(registers[2]);
    edx = static_cast<uint32_t>(registers[3]);
#elif defined(__x86_64__) || defined(__i386__)
    __asm__ volatile(
        "cpuid"
        : "=a"(eax), "=b"(ebx), "=c"(ecx), "=d"(edx)
        : "a"(leaf), "c"(subleaf));
#else
    (void)leaf;
    (void)subleaf;
    eax = ebx = ecx = edx = 0;
#endif
}

static uint64_t Xgetbv0()
{
#if defined(_MSC_VER) && (defined(_M_X64) || defined(_M_IX86))
    return _xgetbv(0);
#elif defined(__x86_64__) || defined(__i386__)
    uint32_t eax, edx;
    __asm__ volatile(
        ".byte 0x0f, 0x01, 0xd0"
        : "=a"(eax), "=d"(edx)
        : "c"(0));
    return (static_cast<uint64_t>(edx) << 32) | eax;
#else
    return 0;
#endif
}

static CpuFeatures DetectCpuFeatures()
{
    CpuFeatures result = { false, false, false, false };
    uint32_t eax, ebx, ecx, edx;
    Cpuid(0, 0, eax, ebx, ecx, edx);
    const uint32_t max_basic = eax;
    Cpuid(0x80000000U, 0, eax, ebx, ecx, edx);
    const uint32_t max_extended = eax;

    if (max_basic < 7)
        return result;

    Cpuid(1, 0, eax, ebx, ecx, edx);
    const bool osxsave = (ecx & (1U << 27)) != 0;
    const bool avx = (ecx & (1U << 28)) != 0;
    if (!osxsave || !avx || (Xgetbv0() & 0xe6U) != 0xe6U)
        return result;

    Cpuid(7, 0, eax, ebx, ecx, edx);
    result.avx2 = (ebx & (1U << 5)) != 0;
    result.avx512f = (ebx & (1U << 16)) != 0;
    result.avx512vl = (ebx & (1U << 31)) != 0;

    if (max_extended >= 0x80000008U)
    {
        Cpuid(0x80000008U, 0, eax, ebx, ecx, edx);
        result.rdpru = (ebx & (1U << 4)) != 0;
    }
    return result;
}

static bool PinToCpu(unsigned requested, unsigned& selected)
{
#if defined(__linux__)
    cpu_set_t allowed;
    CPU_ZERO(&allowed);
    if (sched_getaffinity(0, sizeof(allowed), &allowed) != 0)
        return false;

    if (requested != ~0U)
    {
        if (requested >= CPU_SETSIZE || !CPU_ISSET(requested, &allowed))
            return false;
        selected = requested;
    }
    else
    {
        selected = ~0U;
        for (unsigned cpu = 0; cpu < CPU_SETSIZE; ++cpu)
        {
            if (CPU_ISSET(cpu, &allowed))
            {
                selected = cpu;
                break;
            }
        }
        if (selected == ~0U)
            return false;
    }

    cpu_set_t one;
    CPU_ZERO(&one);
    CPU_SET(selected, &one);
    return sched_setaffinity(0, sizeof(one), &one) == 0;
#else
    (void)requested;
    selected = ~0U;
    return false;
#endif
}

static void PrintUsage(const char* program)
{
    fprintf(stderr, "Usage: %s [--quick] [--verify-only] [--cpu N]\n", program);
}


} // namespace


int main(int argc, char** argv)
{
    bool quick = false;
    bool verify_only = false;
    unsigned requested_cpu = ~0U;
    for (int i = 1; i < argc; ++i)
    {
        if (strcmp(argv[i], "--quick") == 0)
        {
            quick = true;
        }
        else if (strcmp(argv[i], "--verify-only") == 0)
        {
            verify_only = true;
        }
        else if (strcmp(argv[i], "--cpu") == 0 && i + 1 < argc)
        {
            const char* const text = argv[++i];
            char* end = 0;
            errno = 0;
            const unsigned long value = strtoul(text, &end, 10);
            if (errno != 0 || end == text || *end != '\0' ||
                value > 65535UL)
            {
                PrintUsage(argv[0]);
                return 2;
            }
            requested_cpu = static_cast<unsigned>(value);
        }
        else
        {
            PrintUsage(argv[0]);
            return 2;
        }
    }

#if defined(LEO_FF8XOR_XOR3_AVX512_BUILT)
    const bool benchmark_built = true;
#else
    const bool benchmark_built = false;
#endif
    const CpuFeatures features = DetectCpuFeatures();
    printf("FF8 XOR3 / VPTERNLOG microbenchmark\n");
    printf("build_avx512=%s cpu_avx2=%s cpu_avx512f=%s cpu_avx512vl=%s rdpru=%s\n",
        benchmark_built ? "yes" : "no",
        features.avx2 ? "yes" : "no",
        features.avx512f ? "yes" : "no",
        features.avx512vl ? "yes" : "no",
        features.rdpru ? "yes" : "no");

    if (!benchmark_built || !features.avx2 ||
        !features.avx512f || !features.avx512vl)
    {
        printf("SKIP: AVX2 + AVX-512F + AVX-512VL and a built kernel are required.\n");
        return 77;
    }

    unsigned selected_cpu = ~0U;
    if (!PinToCpu(requested_cpu, selected_cpu))
    {
        if (!verify_only)
        {
            fprintf(stderr,
                "Unable to pin benchmark to requested/first allowed CPU.\n");
            return 1;
        }
        printf("pinning=unavailable (allowed for correctness-only run)\n");
    }
    printf("pinned_cpu=%u quick=%s verify_only=%s "
           "measurement=ABBA wall_clock=steady_clock\n",
        selected_cpu, quick ? "yes" : "no", verify_only ? "yes" : "no");
    printf("frequency=RDPRU_APERF_GHz_when_available; TSC_GHz_is_reference_clock\n");

    return ff8xor_xor3_benchmark_run(
        quick, verify_only, selected_cpu, features.rdpru);
}
