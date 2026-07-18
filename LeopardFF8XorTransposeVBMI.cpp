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

#include "LeopardFF8XorTransposeAVX512.h"
#include "LeopardCommon.h"

#if defined(LEO_HAS_FF8) && \
    defined(LEO_FF8XOR_BUILD_AVX512_VBMI_TRANSPOSE)

#include <immintrin.h>

namespace leopard { namespace ff8xor { namespace transpose { namespace avx512 {


#define LEO_PAIR_INDEX(A, B) A, static_cast<uint8_t>(64 + (B))
alignas(64) static const uint8_t kPairLow[64] = {
    LEO_PAIR_INDEX(0,0), LEO_PAIR_INDEX(1,1),
    LEO_PAIR_INDEX(2,2), LEO_PAIR_INDEX(3,3),
    LEO_PAIR_INDEX(4,4), LEO_PAIR_INDEX(5,5),
    LEO_PAIR_INDEX(6,6), LEO_PAIR_INDEX(7,7),
    LEO_PAIR_INDEX(8,8), LEO_PAIR_INDEX(9,9),
    LEO_PAIR_INDEX(10,10), LEO_PAIR_INDEX(11,11),
    LEO_PAIR_INDEX(12,12), LEO_PAIR_INDEX(13,13),
    LEO_PAIR_INDEX(14,14), LEO_PAIR_INDEX(15,15),
    LEO_PAIR_INDEX(16,16), LEO_PAIR_INDEX(17,17),
    LEO_PAIR_INDEX(18,18), LEO_PAIR_INDEX(19,19),
    LEO_PAIR_INDEX(20,20), LEO_PAIR_INDEX(21,21),
    LEO_PAIR_INDEX(22,22), LEO_PAIR_INDEX(23,23),
    LEO_PAIR_INDEX(24,24), LEO_PAIR_INDEX(25,25),
    LEO_PAIR_INDEX(26,26), LEO_PAIR_INDEX(27,27),
    LEO_PAIR_INDEX(28,28), LEO_PAIR_INDEX(29,29),
    LEO_PAIR_INDEX(30,30), LEO_PAIR_INDEX(31,31)
};
alignas(64) static const uint8_t kPairHigh[64] = {
    LEO_PAIR_INDEX(32,32), LEO_PAIR_INDEX(33,33),
    LEO_PAIR_INDEX(34,34), LEO_PAIR_INDEX(35,35),
    LEO_PAIR_INDEX(36,36), LEO_PAIR_INDEX(37,37),
    LEO_PAIR_INDEX(38,38), LEO_PAIR_INDEX(39,39),
    LEO_PAIR_INDEX(40,40), LEO_PAIR_INDEX(41,41),
    LEO_PAIR_INDEX(42,42), LEO_PAIR_INDEX(43,43),
    LEO_PAIR_INDEX(44,44), LEO_PAIR_INDEX(45,45),
    LEO_PAIR_INDEX(46,46), LEO_PAIR_INDEX(47,47),
    LEO_PAIR_INDEX(48,48), LEO_PAIR_INDEX(49,49),
    LEO_PAIR_INDEX(50,50), LEO_PAIR_INDEX(51,51),
    LEO_PAIR_INDEX(52,52), LEO_PAIR_INDEX(53,53),
    LEO_PAIR_INDEX(54,54), LEO_PAIR_INDEX(55,55),
    LEO_PAIR_INDEX(56,56), LEO_PAIR_INDEX(57,57),
    LEO_PAIR_INDEX(58,58), LEO_PAIR_INDEX(59,59),
    LEO_PAIR_INDEX(60,60), LEO_PAIR_INDEX(61,61),
    LEO_PAIR_INDEX(62,62), LEO_PAIR_INDEX(63,63)
};
#undef LEO_PAIR_INDEX

#define LEO_QUAD_INDEX(A) \
    A, A + 1, static_cast<uint8_t>(64 + (A)), \
    static_cast<uint8_t>(65 + (A))
alignas(64) static const uint8_t kQuadLow[64] = {
    LEO_QUAD_INDEX(0), LEO_QUAD_INDEX(2),
    LEO_QUAD_INDEX(4), LEO_QUAD_INDEX(6),
    LEO_QUAD_INDEX(8), LEO_QUAD_INDEX(10),
    LEO_QUAD_INDEX(12), LEO_QUAD_INDEX(14),
    LEO_QUAD_INDEX(16), LEO_QUAD_INDEX(18),
    LEO_QUAD_INDEX(20), LEO_QUAD_INDEX(22),
    LEO_QUAD_INDEX(24), LEO_QUAD_INDEX(26),
    LEO_QUAD_INDEX(28), LEO_QUAD_INDEX(30)
};
alignas(64) static const uint8_t kQuadHigh[64] = {
    LEO_QUAD_INDEX(32), LEO_QUAD_INDEX(34),
    LEO_QUAD_INDEX(36), LEO_QUAD_INDEX(38),
    LEO_QUAD_INDEX(40), LEO_QUAD_INDEX(42),
    LEO_QUAD_INDEX(44), LEO_QUAD_INDEX(46),
    LEO_QUAD_INDEX(48), LEO_QUAD_INDEX(50),
    LEO_QUAD_INDEX(52), LEO_QUAD_INDEX(54),
    LEO_QUAD_INDEX(56), LEO_QUAD_INDEX(58),
    LEO_QUAD_INDEX(60), LEO_QUAD_INDEX(62)
};
#undef LEO_QUAD_INDEX

#define LEO_OCT_INDEX(A) \
    A, A + 1, A + 2, A + 3, static_cast<uint8_t>(64 + (A)), \
    static_cast<uint8_t>(65 + (A)), static_cast<uint8_t>(66 + (A)), \
    static_cast<uint8_t>(67 + (A))
alignas(64) static const uint8_t kOctLow[64] = {
    LEO_OCT_INDEX(0), LEO_OCT_INDEX(4),
    LEO_OCT_INDEX(8), LEO_OCT_INDEX(12),
    LEO_OCT_INDEX(16), LEO_OCT_INDEX(20),
    LEO_OCT_INDEX(24), LEO_OCT_INDEX(28)
};
alignas(64) static const uint8_t kOctHigh[64] = {
    LEO_OCT_INDEX(32), LEO_OCT_INDEX(36),
    LEO_OCT_INDEX(40), LEO_OCT_INDEX(44),
    LEO_OCT_INDEX(48), LEO_OCT_INDEX(52),
    LEO_OCT_INDEX(56), LEO_OCT_INDEX(60)
};
#undef LEO_OCT_INDEX

template <int Shift>
static LEO_FORCE_INLINE void SwapBitGroups(
    __m512i& a, __m512i& b, __m512i mask)
{
    const __m512i bits = _mm512_and_si512(
        _mm512_xor_si512(a, _mm512_srli_epi16(b, Shift)), mask);
    a = _mm512_xor_si512(a, bits);
    b = _mm512_xor_si512(b, _mm512_slli_epi16(bits, Shift));
}

static LEO_FORCE_INLINE void Pair(
    __m512i a, __m512i b, __m512i& low, __m512i& high)
{
    low = _mm512_permutex2var_epi8(
        a, _mm512_load_si512(kPairLow), b);
    high = _mm512_permutex2var_epi8(
        a, _mm512_load_si512(kPairHigh), b);
}

static LEO_FORCE_INLINE void Quad(
    __m512i a, __m512i b, __m512i& low, __m512i& high)
{
    low = _mm512_permutex2var_epi8(
        a, _mm512_load_si512(kQuadLow), b);
    high = _mm512_permutex2var_epi8(
        a, _mm512_load_si512(kQuadHigh), b);
}

static LEO_FORCE_INLINE void StoreOcts(
    __m512i a, __m512i b, uint8_t* low, uint8_t* high)
{
    _mm512_storeu_si512(low,
        _mm512_permutex2var_epi8(
            a, _mm512_load_si512(kOctLow), b));
    _mm512_storeu_si512(high,
        _mm512_permutex2var_epi8(
            a, _mm512_load_si512(kOctHigh), b));
}

bool VbmiKernelsBuilt()
{
    return true;
}

uint64_t PlaneToPackedVbmi(
    const void* plane_void,
    void* packed_void,
    uint64_t plane_bytes,
    uint64_t group)
{
    const uint8_t* plane = static_cast<const uint8_t*>(plane_void);
    uint8_t* packed = static_cast<uint8_t*>(packed_void);

    if (group > plane_bytes)
        return group;

    while (plane_bytes - group >= 64)
    {
        __m512i x0 = _mm512_loadu_si512(
            plane + plane_bytes * 7 + group);
        __m512i x1 = _mm512_loadu_si512(
            plane + plane_bytes * 6 + group);
        __m512i x2 = _mm512_loadu_si512(
            plane + plane_bytes * 5 + group);
        __m512i x3 = _mm512_loadu_si512(
            plane + plane_bytes * 4 + group);
        __m512i x4 = _mm512_loadu_si512(
            plane + plane_bytes * 3 + group);
        __m512i x5 = _mm512_loadu_si512(
            plane + plane_bytes * 2 + group);
        __m512i x6 = _mm512_loadu_si512(
            plane + plane_bytes + group);
        __m512i x7 = _mm512_loadu_si512(plane + group);

        __m512i mask = _mm512_set1_epi8(0x55);
        SwapBitGroups<1>(x0, x1, mask);
        SwapBitGroups<1>(x2, x3, mask);
        SwapBitGroups<1>(x4, x5, mask);
        SwapBitGroups<1>(x6, x7, mask);
        mask = _mm512_set1_epi8(0x33);
        SwapBitGroups<2>(x0, x2, mask);
        SwapBitGroups<2>(x1, x3, mask);
        SwapBitGroups<2>(x4, x6, mask);
        SwapBitGroups<2>(x5, x7, mask);
        mask = _mm512_set1_epi8(0x0f);
        SwapBitGroups<4>(x0, x4, mask);
        SwapBitGroups<4>(x1, x5, mask);
        SwapBitGroups<4>(x2, x6, mask);
        SwapBitGroups<4>(x3, x7, mask);

        __m512i p01a, p01b, p23a, p23b;
        __m512i p45a, p45b, p67a, p67b;
        Pair(x7, x6, p01a, p01b);
        Pair(x5, x4, p23a, p23b);
        Pair(x3, x2, p45a, p45b);
        Pair(x1, x0, p67a, p67b);

        __m512i q03a, q03b, q03c, q03d;
        __m512i q47a, q47b, q47c, q47d;
        Quad(p01a, p23a, q03a, q03b);
        Quad(p01b, p23b, q03c, q03d);
        Quad(p45a, p67a, q47a, q47b);
        Quad(p45b, p67b, q47c, q47d);

        uint8_t* output = packed + group * 8;
        StoreOcts(q03a, q47a, output, output + 64);
        StoreOcts(q03b, q47b, output + 128, output + 192);
        StoreOcts(q03c, q47c, output + 256, output + 320);
        StoreOcts(q03d, q47d, output + 384, output + 448);
        group += 64;
    }
    return group;
}


}}}} // namespace leopard::ff8xor::transpose::avx512

#else

namespace leopard { namespace ff8xor { namespace transpose { namespace avx512 {

bool VbmiKernelsBuilt()
{
    return false;
}

uint64_t PlaneToPackedVbmi(
    const void*, void*, uint64_t, uint64_t start_group)
{
    return start_group;
}

}}}} // namespace leopard::ff8xor::transpose::avx512

#endif
