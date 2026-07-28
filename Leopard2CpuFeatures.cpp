/*
    Copyright (c) 2017 Christopher A. Taylor.  All rights reserved.
    See Leopard2Backend.h for the full BSD license notice.
*/

#include "Leopard2Backend.h"

#include <errno.h>
#include <limits.h>
#include <new>
#include <stdio.h>
#include <string.h>

#if defined(__linux__)
# include <sched.h>
#endif

#if defined(_MSC_VER)
# include <intrin.h>
#endif

namespace leopard { namespace backend {

static const uint32_t kSSSE3 = 0x00000200U;
static const uint32_t kOSXSAVE = 0x08000000U;
static const uint32_t kAVX = 0x10000000U;
static const uint32_t kAVX2 = 0x00000020U;
static const uint32_t kAVX512F = 0x00010000U;
static const uint32_t kAVX512BW = 0x40000000U;
static const uint32_t kAVX512VL = 0x80000000U;
// CPUID.(7,0):ECX bit 8.
static const uint32_t kGFNI = 1U << 8;

#if defined(__linux__)

static const size_t kMaximumAffinityBytes = 16U * 1024U;
static const unsigned kMaximumCacheIndexes = 32;
static const unsigned kMaximumCacheLevel = 16;
static const uint64_t kMinimumPlausibleL3Bytes =
    UINT64_C(1024) * 1024;
static const uint64_t kMaximumPlausibleL3Bytes =
    UINT64_C(1024) * 1024 * 1024 * 1024;

struct LinuxAffinityMask
{
    unsigned long words[
        kMaximumAffinityBytes / sizeof(unsigned long)];
    size_t bytes;
};

enum LinuxTextReadResult
{
    LinuxTextReadSuccess,
    LinuxTextReadMissing,
    LinuxTextReadFailure
};

static bool ParseStrictUnsigned(
    const char* text,
    size_t length,
    uint64_t maximum,
    uint64_t& value)
{
    if (!text || length == 0)
        return false;
    uint64_t parsed = 0;
    for (size_t i = 0; i < length; ++i)
    {
        const unsigned char character =
            static_cast<unsigned char>(text[i]);
        if (character < '0' || character > '9')
            return false;
        const unsigned digit = character - '0';
        if (parsed > (maximum - digit) / 10U)
            return false;
        parsed = parsed * 10U + digit;
    }
    value = parsed;
    return true;
}

static bool StripSingleNewline(const char* text, size_t& length)
{
    if (!text || length == 0)
        return false;
    if (text[length - 1] == '\n')
        --length;
    return length != 0 &&
        (length == 0 || text[length - 1] != '\r') &&
        memchr(text, '\n', length) == NULL &&
        memchr(text, '\r', length) == NULL;
}

static bool ParseLinuxCacheSizeText(
    const char* text,
    uint64_t& bytes)
{
    if (!text)
        return false;
    size_t length = strlen(text);
    if (!StripSingleNewline(text, length))
        return false;

    uint64_t multiplier = 1;
    const char suffix = text[length - 1];
    if (suffix == 'K' || suffix == 'M' || suffix == 'G')
    {
        --length;
        if (length == 0)
            return false;
        multiplier = UINT64_C(1024);
        if (suffix == 'M')
            multiplier *= UINT64_C(1024);
        else if (suffix == 'G')
            multiplier *= UINT64_C(1024) * 1024;
    }

    uint64_t value = 0;
    if (!ParseStrictUnsigned(
            text, length, UINT64_MAX / multiplier, value) ||
        value == 0)
        return false;
    bytes = value * multiplier;
    return bytes <= kMaximumPlausibleL3Bytes;
}

static bool ParseLinuxCacheLevel(
    const char* text,
    unsigned& level)
{
    if (!text)
        return false;
    size_t length = strlen(text);
    if (!StripSingleNewline(text, length))
        return false;
    uint64_t parsed = 0;
    if (!ParseStrictUnsigned(
            text, length, kMaximumCacheLevel, parsed) ||
        parsed == 0)
        return false;
    level = static_cast<unsigned>(parsed);
    return true;
}

static bool ParseLinuxCpuListContains(
    const char* text,
    uint32_t queried_cpu)
{
    if (!text)
        return false;
    size_t length = strlen(text);
    if (!StripSingleNewline(text, length))
        return false;

    size_t cursor = 0;
    uint64_t previous_end = 0;
    bool have_previous = false;
    bool contains = false;
    while (cursor < length)
    {
        const size_t start_begin = cursor;
        while (cursor < length &&
               text[cursor] >= '0' && text[cursor] <= '9')
            ++cursor;
        uint64_t start = 0;
        if (cursor == start_begin ||
            !ParseStrictUnsigned(
                text + start_begin, cursor - start_begin,
                UINT32_MAX, start))
            return false;

        uint64_t end = start;
        if (cursor < length && text[cursor] == '-')
        {
            ++cursor;
            const size_t end_begin = cursor;
            while (cursor < length &&
                   text[cursor] >= '0' && text[cursor] <= '9')
                ++cursor;
            if (cursor == end_begin ||
                !ParseStrictUnsigned(
                    text + end_begin, cursor - end_begin,
                    UINT32_MAX, end) ||
                end < start)
                return false;
        }
        if (have_previous && start <= previous_end)
            return false;
        have_previous = true;
        previous_end = end;
        if (queried_cpu >= start && queried_cpu <= end)
            contains = true;

        if (cursor == length)
            break;
        if (text[cursor] != ',')
            return false;
        ++cursor;
        if (cursor == length)
            return false;
    }
    return contains;
}

static LinuxTextReadResult ReadLinuxTextFile(
    const char* path,
    char* text,
    size_t capacity)
{
    if (!path || !text || capacity < 2)
        return LinuxTextReadFailure;
    errno = 0;
    FILE* file = fopen(path, "rb");
    if (!file)
        return errno == ENOENT || errno == ENOTDIR
            ? LinuxTextReadMissing : LinuxTextReadFailure;

    const size_t length = fread(text, 1, capacity - 1, file);
    bool valid = ferror(file) == 0;
    if (valid && length == capacity - 1)
        valid = fgetc(file) == EOF && feof(file) != 0;
    if (fclose(file) != 0)
        valid = false;
    if (!valid || memchr(text, '\0', length) != NULL)
        return LinuxTextReadFailure;
    text[length] = '\0';
    return LinuxTextReadSuccess;
}

static LinuxTextReadResult ReadLinuxCacheAttribute(
    const char* cpu_root,
    uint32_t cpu,
    unsigned index,
    const char* attribute,
    char* text,
    size_t capacity)
{
    char path[PATH_MAX];
    const int length = snprintf(
        path, sizeof(path), "%s/cpu%u/cache/index%u/%s",
        cpu_root, cpu, index, attribute);
    if (length < 0 || static_cast<size_t>(length) >= sizeof(path))
        return LinuxTextReadFailure;
    return ReadLinuxTextFile(path, text, capacity);
}

static bool LinuxL3BytesForCpu(
    const char* cpu_root,
    uint32_t cpu,
    uint64_t& bytes)
{
    uint64_t smallest_l3_bytes = 0;
    bool found_l3 = false;
    for (unsigned index = 0; index < kMaximumCacheIndexes; ++index)
    {
        char type[32];
        const LinuxTextReadResult type_result = ReadLinuxCacheAttribute(
            cpu_root, cpu, index, "type", type, sizeof(type));
        if (type_result == LinuxTextReadMissing)
        {
            if (index == 0)
                return false;
            break;
        }
        if (type_result != LinuxTextReadSuccess)
            return false;

        size_t type_length = strlen(type);
        if (!StripSingleNewline(type, type_length))
            return false;
        const bool data = type_length == 4 &&
            memcmp(type, "Data", 4) == 0;
        const bool unified = type_length == 7 &&
            memcmp(type, "Unified", 7) == 0;
        const bool instruction = type_length == 11 &&
            memcmp(type, "Instruction", 11) == 0;
        if (!data && !unified && !instruction)
            return false;
        if (instruction)
            continue;

        char level_text[32];
        char size_text[64];
        char shared_cpu_list[512];
        if (ReadLinuxCacheAttribute(
                cpu_root, cpu, index, "level",
                level_text, sizeof(level_text)) != LinuxTextReadSuccess ||
            ReadLinuxCacheAttribute(
                cpu_root, cpu, index, "size",
                size_text, sizeof(size_text)) != LinuxTextReadSuccess ||
            ReadLinuxCacheAttribute(
                cpu_root, cpu, index, "shared_cpu_list",
                shared_cpu_list,
                sizeof(shared_cpu_list)) != LinuxTextReadSuccess)
            return false;

        unsigned level = 0;
        uint64_t cache_bytes = 0;
        if (!ParseLinuxCacheLevel(level_text, level) ||
            !ParseLinuxCacheSizeText(size_text, cache_bytes) ||
            !ParseLinuxCpuListContains(shared_cpu_list, cpu))
            return false;

        if (level == 3 &&
            (!found_l3 || cache_bytes < smallest_l3_bytes))
        {
            found_l3 = true;
            smallest_l3_bytes = cache_bytes;
        }
    }
    if (!found_l3 ||
        smallest_l3_bytes < kMinimumPlausibleL3Bytes ||
        smallest_l3_bytes > kMaximumPlausibleL3Bytes)
        return false;
    bytes = smallest_l3_bytes;
    return true;
}

static bool GetLinuxAffinityMask(LinuxAffinityMask& mask)
{
    memset(&mask, 0, sizeof(mask));
    for (size_t bytes = sizeof(cpu_set_t);
         bytes <= kMaximumAffinityBytes;
         bytes *= 2)
    {
        memset(mask.words, 0, bytes);
        if (sched_getaffinity(
                0, bytes,
                reinterpret_cast<cpu_set_t*>(mask.words)) == 0)
        {
            mask.bytes = bytes;
            return true;
        }
        if (errno != EINVAL || bytes == kMaximumAffinityBytes)
            break;
    }
    return false;
}

static bool SameLinuxAffinityMask(
    const LinuxAffinityMask& left,
    const LinuxAffinityMask& right)
{
    return left.bytes == right.bytes &&
        left.bytes != 0 &&
        memcmp(left.words, right.words, left.bytes) == 0;
}

#ifdef LEO2_ENABLE_TEST_HOOKS
static bool MinimumLinuxL3ForCpuList(
    const char* cpu_root,
    const uint32_t* cpu_ids,
    size_t cpu_count,
    uint64_t& bytes)
{
    if (!cpu_root || !cpu_ids || cpu_count == 0)
        return false;
    uint64_t minimum = UINT64_MAX;
    uint32_t previous = 0;
    for (size_t i = 0; i < cpu_count; ++i)
    {
        if (i != 0 && cpu_ids[i] <= previous)
            return false;
        previous = cpu_ids[i];
        uint64_t cpu_bytes = 0;
        if (!LinuxL3BytesForCpu(cpu_root, cpu_ids[i], cpu_bytes))
            return false;
        if (cpu_bytes < minimum)
            minimum = cpu_bytes;
    }
    if (minimum == UINT64_MAX)
        return false;
    bytes = minimum;
    return true;
}
#endif

static bool MinimumLinuxL3ForAffinity(
    const char* cpu_root,
    const LinuxAffinityMask& affinity,
    uint64_t& bytes)
{
    if (!cpu_root || affinity.bytes == 0)
        return false;
    uint64_t minimum = UINT64_MAX;
    size_t cpu_count = 0;
    const size_t word_count =
        affinity.bytes / sizeof(unsigned long);
    const unsigned bits_per_word =
        static_cast<unsigned>(sizeof(unsigned long) * CHAR_BIT);
    for (size_t word_index = 0;
         word_index < word_count; ++word_index)
    {
        const unsigned long word = affinity.words[word_index];
        for (unsigned bit = 0; bit < bits_per_word; ++bit)
        {
            if ((word & (static_cast<unsigned long>(1) << bit)) == 0)
                continue;
            const uint64_t cpu_index =
                static_cast<uint64_t>(word_index) * bits_per_word + bit;
            if (cpu_index > UINT32_MAX)
                return false;
            uint64_t cpu_bytes = 0;
            if (!LinuxL3BytesForCpu(
                    cpu_root, static_cast<uint32_t>(cpu_index),
                    cpu_bytes))
                return false;
            if (cpu_bytes < minimum)
                minimum = cpu_bytes;
            ++cpu_count;
        }
    }
    if (cpu_count == 0 || minimum == UINT64_MAX)
        return false;
    bytes = minimum;
    return true;
}

#endif // __linux__

X86ProcessorIdentity ClassifyX86Processor(
    uint32_t vendor_ebx,
    uint32_t vendor_edx,
    uint32_t vendor_ecx,
    uint32_t leaf1_eax)
{
    static const uint32_t kAuthenticAMDEbx = 0x68747541U;
    static const uint32_t kAuthenticAMDEdx = 0x69746e65U;
    static const uint32_t kAuthenticAMDEcx = 0x444d4163U;
    X86ProcessorIdentity identity = { false, 0, 0 };
    identity.authentic_amd =
        vendor_ebx == kAuthenticAMDEbx &&
        vendor_edx == kAuthenticAMDEdx &&
        vendor_ecx == kAuthenticAMDEcx;

    const uint32_t base_family = (leaf1_eax >> 8) & 0xfU;
    const uint32_t extended_family = (leaf1_eax >> 20) & 0xffU;
    const uint32_t base_model = (leaf1_eax >> 4) & 0xfU;
    const uint32_t extended_model = (leaf1_eax >> 16) & 0xfU;
    identity.family = static_cast<uint16_t>(base_family == 0xfU
        ? base_family + extended_family : base_family);
    identity.model = static_cast<uint16_t>(
        base_family == 0x6U || base_family == 0xfU
            ? base_model + (extended_model << 4) : base_model);
    return identity;
}

X86Features ClassifyX86Features(
    uint32_t maximum_basic_leaf,
    uint32_t leaf1_ecx,
    uint32_t leaf7_ebx,
    uint32_t leaf7_ecx,
    uint64_t xcr0)
{
    X86Features features = { false, false, false, false };
    if (maximum_basic_leaf < 1)
        return features;

    features.ssse3 = (leaf1_ecx & kSSSE3) != 0;
    const uint32_t avx_os_mask = kAVX | kOSXSAVE;
    const bool ymm_enabled = (leaf1_ecx & avx_os_mask) == avx_os_mask &&
        (xcr0 & 0x6U) == 0x6U;
    features.avx2 = maximum_basic_leaf >= 7 && ymm_enabled &&
        (leaf7_ebx & kAVX2) != 0;
    const uint32_t avx512_mask = kAVX512F | kAVX512BW | kAVX512VL;
    const bool zmm_enabled = (xcr0 & 0xe6U) == 0xe6U;
    features.avx512 = features.avx2 && zmm_enabled &&
        (leaf7_ebx & avx512_mask) == avx512_mask;
    // GFNI is CPUID.(7,0):ECX bit 8.  The VEX-encoded 256-bit forms the GFNI
    // backend emits require exactly GFNI plus AVX2-class YMM state, so the
    // qualification is the AVX2 gate with the additional feature bit; ZMM
    // state is deliberately not required.
    features.gfni = features.avx2 && (leaf7_ecx & kGFNI) != 0;
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
        X86Features features = { false, false, false, false };
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
    uint32_t leaf7_ecx = 0;
    if (maximum_basic_leaf >= 7)
    {
        ReadCPUID(7, registers);
        leaf7_ebx = registers[1];
        leaf7_ecx = registers[2];
    }
    return ClassifyX86Features(
        maximum_basic_leaf, leaf1_ecx, leaf7_ebx, leaf7_ecx, xcr0);
}

X86ProcessorIdentity DetectX86Processor()
{
    if (!CPUIDSupported())
    {
        X86ProcessorIdentity identity = { false, 0, 0 };
        return identity;
    }
    uint32_t registers[4] = {};
    ReadCPUID(0, registers);
    const uint32_t maximum_basic_leaf = registers[0];
    const uint32_t vendor_ebx = registers[1];
    const uint32_t vendor_ecx = registers[2];
    const uint32_t vendor_edx = registers[3];
    uint32_t leaf1_eax = 0;
    if (maximum_basic_leaf >= 1)
    {
        ReadCPUID(1, registers);
        leaf1_eax = registers[0];
    }
    return ClassifyX86Processor(
        vendor_ebx, vendor_edx, vendor_ecx, leaf1_eax);
}

#else

X86Features DetectX86Features()
{
    X86Features features = { false, false, false, false };
    return features;
}

X86ProcessorIdentity DetectX86Processor()
{
    X86ProcessorIdentity identity = { false, 0, 0 };
    return identity;
}

#endif

bool IsCalibratedAutoAVX512EncodeProcessor(
    const X86ProcessorIdentity& identity)
{
    // Family 1Ah, model 44h is the Zen 5 Granite Ridge class on which the
    // complete-codec table was calibrated across both cache complexes.
    return identity.authentic_amd && identity.family == 0x1aU &&
        identity.model == 0x44U;
}

bool IsCalibratedAutoAVX512EncodeHost()
{
    return IsCalibratedAutoAVX512EncodeProcessor(DetectX86Processor());
}

bool IsCalibratedAutoGFNIProcessor(const X86ProcessorIdentity& identity)
{
    /*
        Same Zen 5 Granite Ridge pin as the AVX-512 encode exception.  The
        qualifying evidence is the 2026-07-28 44-cell same-binary screen
        (both fields, K=2..4096, low and high rate, 1 KiB-1 MiB): zero digest
        mismatches, zero encode regressions, and the single flagged decode
        cell re-measured at parity over four clean trials.  This predicate
        gates a default-off candidate only; flipping the AUTO default
        additionally requires the isolated exact-main campaign
        (docs/leopard2_gfni_codec.md requirement 4).
    */
    return identity.authentic_amd && identity.family == 0x1aU &&
        identity.model == 0x44U;
}

bool IsCalibratedAutoGFNIHost()
{
    return IsCalibratedAutoGFNIProcessor(DetectX86Processor());
}

uint64_t DetectConservativeL3Bytes()
{
#if defined(__linux__)
    static const char kLinuxCpuRoot[] = "/sys/devices/system/cpu";
    LinuxAffinityMask* masks =
        new (std::nothrow) LinuxAffinityMask[2];
    if (!masks)
        return 0;
    uint64_t detected = 0;
    /*
        Retry once if a concurrent affinity change is observed.  A stable
        before/after mask is required because silently combining topology from
        two masks could choose the larger cache on an asymmetric processor.
    */
    for (unsigned attempt = 0; attempt < 2; ++attempt)
    {
        uint64_t minimum = 0;
        if (!GetLinuxAffinityMask(masks[0]) ||
            !MinimumLinuxL3ForAffinity(
                kLinuxCpuRoot, masks[0], minimum) ||
            !GetLinuxAffinityMask(masks[1]))
            break;
        if (SameLinuxAffinityMask(masks[0], masks[1]))
        {
            detected = minimum;
            break;
        }
    }
    delete[] masks;
    return detected;
#endif
    return 0;
}

#ifdef LEO2_ENABLE_TEST_HOOKS
bool TestOnlyParseLinuxCacheSize(
    const char* text,
    uint64_t* bytes_out)
{
#if defined(__linux__)
    if (!bytes_out)
        return false;
    uint64_t bytes = 0;
    if (!ParseLinuxCacheSizeText(text, bytes))
        return false;
    *bytes_out = bytes;
    return true;
#else
    (void)text;
    (void)bytes_out;
    return false;
#endif
}

bool TestOnlyDetectLinuxL3Bytes(
    const char* cpu_root,
    const uint32_t* cpu_ids,
    size_t cpu_count,
    uint64_t* bytes_out)
{
#if defined(__linux__)
    if (!bytes_out)
        return false;
    uint64_t bytes = 0;
    if (!MinimumLinuxL3ForCpuList(
            cpu_root, cpu_ids, cpu_count, bytes))
        return false;
    *bytes_out = bytes;
    return true;
#else
    (void)cpu_root;
    (void)cpu_ids;
    (void)cpu_count;
    (void)bytes_out;
    return false;
#endif
}
#endif

}} // namespace leopard::backend
