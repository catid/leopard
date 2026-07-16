/*
    Copyright (c) 2017 Christopher A. Taylor.  All rights reserved.
    See Leopard2Backend.h for the full BSD license notice.
*/

#include "Leopard2Backend.h"

#if defined(_MSC_VER)
# include <intrin.h>
#endif

namespace leopard { namespace backend {

static const uint32_t kSSSE3 = 0x00000200U;
static const uint32_t kOSXSAVE = 0x08000000U;
static const uint32_t kAVX = 0x10000000U;
static const uint32_t kAVX2 = 0x00000020U;

X86Features ClassifyX86Features(
    uint32_t maximum_basic_leaf,
    uint32_t leaf1_ecx,
    uint32_t leaf7_ebx,
    uint64_t xcr0)
{
    X86Features features = { false, false };
    if (maximum_basic_leaf < 1)
        return features;

    features.ssse3 = (leaf1_ecx & kSSSE3) != 0;
    const uint32_t avx_os_mask = kAVX | kOSXSAVE;
    const bool ymm_enabled = (leaf1_ecx & avx_os_mask) == avx_os_mask &&
        (xcr0 & 0x6U) == 0x6U;
    features.avx2 = maximum_basic_leaf >= 7 && ymm_enabled &&
        (leaf7_ebx & kAVX2) != 0;
    return features;
}

#if defined(__i386__) || defined(__x86_64__) || \
    defined(_M_IX86) || defined(_M_X64) || defined(_M_AMD64)

#if defined(__i386__)
static bool CPUIDSupported()
{
    uint32_t changed = 0;
    uint32_t original = 0;
    __asm__ __volatile__(
        "pushfl; pushfl; popl %0; movl %0, %1; xorl %2, %0; "
        "pushl %0; popfl; pushfl; popl %0; popfl"
        : "=&r"(changed), "=&r"(original)
        : "i"(0x200000));
    return ((changed ^ original) & 0x200000U) != 0;
}
#else
static bool CPUIDSupported()
{
    // CPUID is part of the x86-64 architecture and Windows 7's supported x86
    // processor floor.  The explicit EFLAGS probe above preserves the older
    // GCC i386 path without adding a post-baseline instruction.
    return true;
}
#endif

static void ReadCPUID(uint32_t leaf, uint32_t registers[4])
{
#if defined(_MSC_VER)
    int values[4] = {};
    __cpuidex(values, static_cast<int>(leaf), 0);
    for (unsigned i = 0; i < 4; ++i)
        registers[i] = static_cast<uint32_t>(values[i]);
#elif defined(__i386__) && defined(__PIC__)
    __asm__ __volatile__(
        "xchgl %%ebx, %1; cpuid; xchgl %%ebx, %1"
        : "=a"(registers[0]), "=&r"(registers[1]),
          "=c"(registers[2]), "=d"(registers[3])
        : "0"(leaf), "2"(0U));
#elif defined(__x86_64__) && defined(__PIC__)
    __asm__ __volatile__(
        "xchgq %%rbx, %q1; cpuid; xchgq %%rbx, %q1"
        : "=a"(registers[0]), "=&r"(registers[1]),
          "=c"(registers[2]), "=d"(registers[3])
        : "0"(leaf), "2"(0U));
#else
    __asm__ __volatile__(
        "cpuid"
        : "=a"(registers[0]), "=b"(registers[1]),
          "=c"(registers[2]), "=d"(registers[3])
        : "0"(leaf), "2"(0U));
#endif
}

static uint64_t ReadXCR0()
{
#if defined(_MSC_VER)
    return static_cast<uint64_t>(_xgetbv(0));
#else
    uint32_t eax = 0;
    uint32_t edx = 0;
    __asm__ __volatile__("xgetbv" : "=a"(eax), "=d"(edx) : "c"(0));
    return (static_cast<uint64_t>(edx) << 32) | eax;
#endif
}

X86Features DetectX86Features()
{
    if (!CPUIDSupported())
    {
        X86Features features = { false, false };
        return features;
    }
    uint32_t registers[4] = {};
    ReadCPUID(0, registers);
    const uint32_t maximum_basic_leaf = registers[0];
    uint32_t leaf1_ecx = 0;
    uint32_t leaf7_ebx = 0;
    uint64_t xcr0 = 0;
    if (maximum_basic_leaf >= 1)
    {
        ReadCPUID(1, registers);
        leaf1_ecx = registers[2];
        if ((leaf1_ecx & (kAVX | kOSXSAVE)) == (kAVX | kOSXSAVE))
            xcr0 = ReadXCR0();
    }
    if (maximum_basic_leaf >= 7)
    {
        ReadCPUID(7, registers);
        leaf7_ebx = registers[1];
    }
    return ClassifyX86Features(
        maximum_basic_leaf, leaf1_ecx, leaf7_ebx, xcr0);
}

#else

X86Features DetectX86Features()
{
    X86Features features = { false, false };
    return features;
}

#endif

}} // namespace leopard::backend
