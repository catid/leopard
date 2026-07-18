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
#include "LeopardFF8XorTransposeAVX2.h"

#if defined(LEO_HAS_FF8) && defined(LEO_FF8XOR_BUILD_AVX2_TRANSPOSE)

#include <immintrin.h>

namespace leopard { namespace ff8xor { namespace transpose { namespace avx2 {


static const uint64_t kGroupsPerBlock = 32;
static const uint64_t kPackedBytesPerVector = 32;

static LEO_FORCE_INLINE __m256i Load(const uint8_t* source)
{
    return _mm256_loadu_si256(reinterpret_cast<const __m256i*>(source));
}

static LEO_FORCE_INLINE void Store(uint8_t* destination, __m256i value)
{
    _mm256_storeu_si256(reinterpret_cast<__m256i*>(destination), value);
}

template <int Shift>
static LEO_FORCE_INLINE int GatherHighBits(__m256i value)
{
    // The shift is at most seven.  Within each 16-bit lane, bits shifted out
    // of the low byte cannot reach the high byte's sign bit, so movemask
    // independently extracts the selected bit from all 32 input bytes.
    return _mm256_movemask_epi8(_mm256_slli_epi16(value, Shift));
}

template <int Shift>
static LEO_FORCE_INLINE __m256i GatherPlane(
    __m256i packed0,
    __m256i packed1,
    __m256i packed2,
    __m256i packed3,
    __m256i packed4,
    __m256i packed5,
    __m256i packed6,
    __m256i packed7)
{
    return _mm256_setr_epi32(
        GatherHighBits<Shift>(packed0),
        GatherHighBits<Shift>(packed1),
        GatherHighBits<Shift>(packed2),
        GatherHighBits<Shift>(packed3),
        GatherHighBits<Shift>(packed4),
        GatherHighBits<Shift>(packed5),
        GatherHighBits<Shift>(packed6),
        GatherHighBits<Shift>(packed7));
}

static LEO_FORCE_INLINE void PackedToPlaneBlock(
    const uint8_t* packed,
    uint8_t* plane,
    uint64_t plane_bytes,
    uint64_t group)
{
    const uint8_t* block = packed + group * 8;
    const __m256i packed0 = Load(block);
    const __m256i packed1 = Load(block + kPackedBytesPerVector);
    const __m256i packed2 = Load(block + kPackedBytesPerVector * 2);
    const __m256i packed3 = Load(block + kPackedBytesPerVector * 3);
    const __m256i packed4 = Load(block + kPackedBytesPerVector * 4);
    const __m256i packed5 = Load(block + kPackedBytesPerVector * 5);
    const __m256i packed6 = Load(block + kPackedBytesPerVector * 6);
    const __m256i packed7 = Load(block + kPackedBytesPerVector * 7);

    Store(plane + group,
        GatherPlane<7>(packed0, packed1, packed2, packed3,
            packed4, packed5, packed6, packed7));
    Store(plane + plane_bytes + group,
        GatherPlane<6>(packed0, packed1, packed2, packed3,
            packed4, packed5, packed6, packed7));
    Store(plane + plane_bytes * 2 + group,
        GatherPlane<5>(packed0, packed1, packed2, packed3,
            packed4, packed5, packed6, packed7));
    Store(plane + plane_bytes * 3 + group,
        GatherPlane<4>(packed0, packed1, packed2, packed3,
            packed4, packed5, packed6, packed7));
    Store(plane + plane_bytes * 4 + group,
        GatherPlane<3>(packed0, packed1, packed2, packed3,
            packed4, packed5, packed6, packed7));
    Store(plane + plane_bytes * 5 + group,
        GatherPlane<2>(packed0, packed1, packed2, packed3,
            packed4, packed5, packed6, packed7));
    Store(plane + plane_bytes * 6 + group,
        GatherPlane<1>(packed0, packed1, packed2, packed3,
            packed4, packed5, packed6, packed7));
    Store(plane + plane_bytes * 7 + group,
        GatherPlane<0>(packed0, packed1, packed2, packed3,
            packed4, packed5, packed6, packed7));
}

template <int Shift>
static LEO_FORCE_INLINE void SwapBitGroups(
    __m256i& a,
    __m256i& b,
    __m256i mask)
{
    const __m256i bits = _mm256_and_si256(
        _mm256_xor_si256(a, _mm256_srli_epi16(b, Shift)), mask);
    a = _mm256_xor_si256(a, bits);
    b = _mm256_xor_si256(b, _mm256_slli_epi16(bits, Shift));
}

static LEO_FORCE_INLINE void UnpackBytes(__m256i& a, __m256i& b)
{
    const __m256i old_a = a;
    a = _mm256_unpacklo_epi8(a, b);
    b = _mm256_unpackhi_epi8(old_a, b);
}

static LEO_FORCE_INLINE void UnpackWords(__m256i& a, __m256i& b)
{
    const __m256i old_a = a;
    a = _mm256_unpacklo_epi16(a, b);
    b = _mm256_unpackhi_epi16(old_a, b);
}

static LEO_FORCE_INLINE void UnpackDoubleWords(__m256i& a, __m256i& b)
{
    const __m256i old_a = a;
    a = _mm256_unpacklo_epi32(a, b);
    b = _mm256_unpackhi_epi32(old_a, b);
}

static LEO_FORCE_INLINE void PlaneToPackedBlock(
    const uint8_t* plane,
    uint8_t* packed,
    uint64_t plane_bytes,
    uint64_t group)
{
    // Reverse the row order on input and output around this standard 8x8 bit
    // transpose network.  The result is laneN[group] bit J =
    // planeJ[group] bit N for 32 independent groups at once.
    __m256i x0 = Load(plane + plane_bytes * 7 + group);
    __m256i x1 = Load(plane + plane_bytes * 6 + group);
    __m256i x2 = Load(plane + plane_bytes * 5 + group);
    __m256i x3 = Load(plane + plane_bytes * 4 + group);
    __m256i x4 = Load(plane + plane_bytes * 3 + group);
    __m256i x5 = Load(plane + plane_bytes * 2 + group);
    __m256i x6 = Load(plane + plane_bytes + group);
    __m256i x7 = Load(plane + group);

    __m256i mask = _mm256_set1_epi8(0x55);
    SwapBitGroups<1>(x0, x1, mask);
    SwapBitGroups<1>(x2, x3, mask);
    SwapBitGroups<1>(x4, x5, mask);
    SwapBitGroups<1>(x6, x7, mask);
    mask = _mm256_set1_epi8(0x33);
    SwapBitGroups<2>(x0, x2, mask);
    SwapBitGroups<2>(x1, x3, mask);
    SwapBitGroups<2>(x4, x6, mask);
    SwapBitGroups<2>(x5, x7, mask);
    mask = _mm256_set1_epi8(0x0f);
    SwapBitGroups<4>(x0, x4, mask);
    SwapBitGroups<4>(x1, x5, mask);
    SwapBitGroups<4>(x2, x6, mask);
    SwapBitGroups<4>(x3, x7, mask);

    // Reverse the vector order without adding live aliases, then interleave
    // 8 rows x 32 group bytes in place.  Each stage needs just one temporary;
    // this avoids spilling the eight row registers on AVX2's 16-register ISA.
    __m256i temporary = x0;
    x0 = x7;
    x7 = temporary;
    temporary = x1;
    x1 = x6;
    x6 = temporary;
    temporary = x2;
    x2 = x5;
    x5 = temporary;
    temporary = x3;
    x3 = x4;
    x4 = temporary;

    UnpackBytes(x0, x1);
    UnpackBytes(x2, x3);
    UnpackBytes(x4, x5);
    UnpackBytes(x6, x7);
    UnpackWords(x0, x2);
    UnpackWords(x1, x3);
    UnpackWords(x4, x6);
    UnpackWords(x5, x7);
    UnpackDoubleWords(x0, x4);
    UnpackDoubleWords(x2, x6);
    UnpackDoubleWords(x1, x5);
    UnpackDoubleWords(x3, x7);

    uint8_t* block = packed + group * 8;
    Store(block,
        _mm256_permute2x128_si256(x0, x4, 0x20));
    Store(block + kPackedBytesPerVector,
        _mm256_permute2x128_si256(x2, x6, 0x20));
    Store(block + kPackedBytesPerVector * 2,
        _mm256_permute2x128_si256(x1, x5, 0x20));
    Store(block + kPackedBytesPerVector * 3,
        _mm256_permute2x128_si256(x3, x7, 0x20));
    Store(block + kPackedBytesPerVector * 4,
        _mm256_permute2x128_si256(x0, x4, 0x31));
    Store(block + kPackedBytesPerVector * 5,
        _mm256_permute2x128_si256(x2, x6, 0x31));
    Store(block + kPackedBytesPerVector * 6,
        _mm256_permute2x128_si256(x1, x5, 0x31));
    Store(block + kPackedBytesPerVector * 7,
        _mm256_permute2x128_si256(x3, x7, 0x31));
}

bool KernelsBuilt()
{
    return true;
}

uint64_t PackedToPlane(
    const void* packed_void,
    void* plane_void,
    uint64_t plane_bytes,
    uint64_t group)
{
    const uint8_t* packed = static_cast<const uint8_t*>(packed_void);
    uint8_t* plane = static_cast<uint8_t*>(plane_void);

    if (group > plane_bytes)
        return group;

    while (plane_bytes - group >= kGroupsPerBlock)
    {
        PackedToPlaneBlock(packed, plane, plane_bytes, group);
        group += kGroupsPerBlock;
    }
    return group;
}

uint64_t PlaneToPacked(
    const void* plane_void,
    void* packed_void,
    uint64_t plane_bytes,
    uint64_t group)
{
    const uint8_t* plane = static_cast<const uint8_t*>(plane_void);
    uint8_t* packed = static_cast<uint8_t*>(packed_void);

    if (group > plane_bytes)
        return group;

    while (plane_bytes - group >= kGroupsPerBlock)
    {
        PlaneToPackedBlock(plane, packed, plane_bytes, group);
        group += kGroupsPerBlock;
    }
    return group;
}


}}}} // namespace leopard::ff8xor::transpose::avx2

#else

namespace leopard { namespace ff8xor { namespace transpose { namespace avx2 {


bool KernelsBuilt()
{
    return false;
}

uint64_t PackedToPlane(
    const void*,
    void*,
    uint64_t,
    uint64_t start_group)
{
    return start_group;
}

uint64_t PlaneToPacked(
    const void*,
    void*,
    uint64_t,
    uint64_t start_group)
{
    return start_group;
}


}}}} // namespace leopard::ff8xor::transpose::avx2

#endif
