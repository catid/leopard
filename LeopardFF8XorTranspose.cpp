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

#include "LeopardFF8XorTranspose.h"

#include "LeopardFF8Xor.h"
#include "LeopardFF8XorTransposeAVX2.h"
#include "LeopardFF8XorTransposeAVX512.h"

#include <stddef.h>
#include <string.h>

#if defined(_MSC_VER) && \
    (defined(_M_X64) || defined(_M_AMD64) || defined(_M_IX86))
    #include <intrin.h>
#elif defined(__i386__) || defined(__x86_64__)
    #include <cpuid.h>
#endif

namespace leopard { namespace ff8xor { namespace transpose {


static uint64_t TransposeWord8x8(uint64_t value)
{
    uint64_t temp = (value ^ (value >> 7)) &
        UINT64_C(0x00aa00aa00aa00aa);
    value ^= temp ^ (temp << 7);
    temp = (value ^ (value >> 14)) & UINT64_C(0x0000cccc0000cccc);
    value ^= temp ^ (temp << 14);
    temp = (value ^ (value >> 28)) & UINT64_C(0x00000000f0f0f0f0);
    return value ^ temp ^ (temp << 28);
}

static uint64_t LoadPackedWord(const uint8_t* source)
{
    uint64_t value;
    memcpy(&value, source, sizeof(value));
#if defined(__BYTE_ORDER__) && __BYTE_ORDER__ == __ORDER_BIG_ENDIAN__
    value = __builtin_bswap64(value);
#endif
    return value;
}

static void StorePackedWord(uint8_t* destination, uint64_t value)
{
#if defined(__BYTE_ORDER__) && __BYTE_ORDER__ == __ORDER_BIG_ENDIAN__
    value = __builtin_bswap64(value);
#endif
    memcpy(destination, &value, sizeof(value));
}

static void PackedToPlanePortableTail(
    const uint8_t* packed,
    uint8_t* plane,
    uint64_t plane_bytes,
    uint64_t group)
{
    for (; group < plane_bytes; ++group)
    {
        const uint64_t word = TransposeWord8x8(
            LoadPackedWord(packed + group * 8));
        for (unsigned bit = 0; bit < 8; ++bit)
        {
            plane[bit * plane_bytes + group] =
                static_cast<uint8_t>(word >> (bit * 8));
        }
    }
}

static bool HasAVX512State(
    uint32_t maximum_basic_leaf,
    uint32_t leaf1_ecx,
    uint32_t leaf7_ebx,
    uint64_t xcr0)
{
    const uint32_t required_leaf1 =
        (UINT32_C(1) << 27) | // OSXSAVE
        (UINT32_C(1) << 28);  // AVX
    const uint32_t required_leaf7 =
        (UINT32_C(1) << 16) | // AVX-512F
        (UINT32_C(1) << 30);  // AVX-512BW
    return maximum_basic_leaf >= 7 &&
        (leaf1_ecx & required_leaf1) == required_leaf1 &&
        (leaf7_ebx & required_leaf7) == required_leaf7 &&
        (xcr0 & UINT64_C(0xe6)) == UINT64_C(0xe6);
}

bool IsAVX512BitalgSupported(
    uint32_t maximum_basic_leaf,
    uint32_t leaf1_ecx,
    uint32_t leaf7_ebx,
    uint32_t leaf7_ecx,
    uint64_t xcr0)
{
    return HasAVX512State(maximum_basic_leaf, leaf1_ecx, leaf7_ebx, xcr0) &&
        (leaf7_ecx & (UINT32_C(1) << 12)) != 0; // AVX-512BITALG
}

bool IsAVX512VbmiSupported(
    uint32_t maximum_basic_leaf,
    uint32_t leaf1_ecx,
    uint32_t leaf7_ebx,
    uint32_t leaf7_ecx,
    uint64_t xcr0)
{
    return HasAVX512State(maximum_basic_leaf, leaf1_ecx, leaf7_ebx, xcr0) &&
        (leaf7_ecx & (UINT32_C(1) << 1)) != 0; // AVX-512VBMI
}

struct AVX512TransposeCapabilities
{
    bool Bitalg;
    bool Vbmi;
};

#if defined(_MSC_VER)
    #define LEO_FF8XOR_TRANSPOSE_NOINLINE __declspec(noinline)
#elif defined(__GNUC__) || defined(__clang__)
    #define LEO_FF8XOR_TRANSPOSE_NOINLINE __attribute__((noinline))
#else
    #define LEO_FF8XOR_TRANSPOSE_NOINLINE
#endif

static LEO_FF8XOR_TRANSPOSE_NOINLINE AVX512TransposeCapabilities
DetectAVX512TransposeCapabilities()
{
    AVX512TransposeCapabilities capabilities = { false, false };

#if (defined(LEO_FF8XOR_HAS_AVX512_BITALG_TRANSPOSE) || \
     defined(LEO_FF8XOR_HAS_AVX512_VBMI_TRANSPOSE)) && \
    defined(_MSC_VER) && \
    (defined(_M_X64) || defined(_M_AMD64) || defined(_M_IX86))
    int registers[4] = { 0, 0, 0, 0 };
    __cpuidex(registers, 0, 0);
    const uint32_t maximum_basic_leaf =
        static_cast<uint32_t>(registers[0]);
    if (maximum_basic_leaf < 7)
        return capabilities;
    __cpuidex(registers, 1, 0);
    const uint32_t leaf1_ecx = static_cast<uint32_t>(registers[2]);
    const uint32_t required_leaf1 =
        (UINT32_C(1) << 27) | (UINT32_C(1) << 28);
    if ((leaf1_ecx & required_leaf1) != required_leaf1)
        return capabilities;
    const uint64_t xcr0 = static_cast<uint64_t>(_xgetbv(0));
    __cpuidex(registers, 7, 0);
    const uint32_t leaf7_ebx = static_cast<uint32_t>(registers[1]);
    const uint32_t leaf7_ecx = static_cast<uint32_t>(registers[2]);
#elif (defined(LEO_FF8XOR_HAS_AVX512_BITALG_TRANSPOSE) || \
       defined(LEO_FF8XOR_HAS_AVX512_VBMI_TRANSPOSE)) && \
      (defined(__i386__) || defined(__x86_64__))
    const uint32_t maximum_basic_leaf = __get_cpuid_max(0, NULL);
    if (maximum_basic_leaf < 7)
        return capabilities;
    unsigned eax = 0, ebx = 0, ecx = 0, edx = 0;
    __cpuid_count(1, 0, eax, ebx, ecx, edx);
    const uint32_t leaf1_ecx = ecx;
    const uint32_t required_leaf1 =
        (UINT32_C(1) << 27) | (UINT32_C(1) << 28);
    if ((leaf1_ecx & required_leaf1) != required_leaf1)
        return capabilities;
    uint32_t xcr0_low = 0, xcr0_high = 0;
    __asm__ __volatile__(
        ".byte 0x0f, 0x01, 0xd0"
        : "=a"(xcr0_low), "=d"(xcr0_high) : "c"(0));
    const uint64_t xcr0 = static_cast<uint64_t>(xcr0_low) |
        (static_cast<uint64_t>(xcr0_high) << 32);
    __cpuid_count(7, 0, eax, ebx, ecx, edx);
    const uint32_t leaf7_ebx = ebx;
    const uint32_t leaf7_ecx = ecx;
#else
    return capabilities;
#endif

#ifdef LEO_FF8XOR_HAS_AVX512_BITALG_TRANSPOSE
    capabilities.Bitalg = IsAVX512BitalgSupported(
        maximum_basic_leaf, leaf1_ecx, leaf7_ebx, leaf7_ecx, xcr0) &&
        avx512::BitalgKernelsBuilt();
#endif
#ifdef LEO_FF8XOR_HAS_AVX512_VBMI_TRANSPOSE
    capabilities.Vbmi = IsAVX512VbmiSupported(
        maximum_basic_leaf, leaf1_ecx, leaf7_ebx, leaf7_ecx, xcr0) &&
        avx512::VbmiKernelsBuilt();
#endif
    return capabilities;
}

#undef LEO_FF8XOR_TRANSPOSE_NOINLINE

static const AVX512TransposeCapabilities&
GetAVX512TransposeCapabilities()
{
    static const AVX512TransposeCapabilities capabilities =
        DetectAVX512TransposeCapabilities();
    return capabilities;
}

bool IsModeAvailable(Mode mode)
{
    if (mode == Mode::Avx2)
    {
        return ff8xor::IsKernelModeAvailable(KernelMode::Avx2) &&
            avx2::KernelsBuilt();
    }
    if (mode == Mode::Avx512Bitalg)
        return GetAVX512TransposeCapabilities().Bitalg;
    if (mode == Mode::Avx512Vbmi)
        return GetAVX512TransposeCapabilities().Vbmi;
    return true;
}

void PackedToPlane(
    const void* packed_void,
    void* plane_void,
    uint64_t buffer_bytes,
    Mode mode)
{
    const uint8_t* packed = static_cast<const uint8_t*>(packed_void);
    uint8_t* plane = static_cast<uint8_t*>(plane_void);
    const uint64_t plane_bytes = buffer_bytes / 8;
    uint64_t group = 0;

    if (plane_bytes >= 8 &&
        (mode == Mode::Auto || mode == Mode::Avx512Bitalg) &&
        IsModeAvailable(Mode::Avx512Bitalg))
    {
        group = avx512::PackedToPlaneBitalg(
            packed, plane, plane_bytes, group);
    }

    // BITALG consumes every complete eight-group block.  If it advanced the
    // cursor, fewer than eight groups remain and the 32-group AVX2 kernel can
    // no longer do useful work.  Keep AVX2 dispatch for Auto without BITALG
    // and for an explicitly forced AVX2 control.
    if (group == 0 && plane_bytes >= 32 &&
        (mode == Mode::Auto || mode == Mode::Avx2) &&
        IsModeAvailable(Mode::Avx2))
    {
        group = avx2::PackedToPlane(
            packed, plane, plane_bytes, group);
    }

    PackedToPlanePortableTail(
        packed, plane, plane_bytes, group);
}

void PlaneToPacked(
    const void* plane_void,
    void* packed_void,
    uint64_t buffer_bytes,
    Mode mode)
{
    const uint8_t* plane = static_cast<const uint8_t*>(plane_void);
    uint8_t* packed = static_cast<uint8_t*>(packed_void);
    const uint64_t plane_bytes = buffer_bytes / 8;
    uint64_t group = 0;

    if (plane_bytes >= 64 &&
        (mode == Mode::Auto || mode == Mode::Avx512Vbmi) &&
        IsModeAvailable(Mode::Avx512Vbmi))
    {
        group = avx512::PlaneToPackedVbmi(
            plane, packed, plane_bytes, group);
    }

    if (group <= plane_bytes && plane_bytes - group >= 32 &&
        (mode == Mode::Auto || mode == Mode::Avx2) &&
        IsModeAvailable(Mode::Avx2))
    {
        group = avx2::PlaneToPacked(
            plane, packed, plane_bytes, group);
    }

    for (; group < plane_bytes; ++group)
    {
        uint64_t word = 0;
        for (unsigned bit = 0; bit < 8; ++bit)
        {
            word |= static_cast<uint64_t>(
                plane[bit * plane_bytes + group]) << (bit * 8);
        }
        StorePackedWord(
            packed + group * 8, TransposeWord8x8(word));
    }
}


}}} // namespace leopard::ff8xor::transpose
