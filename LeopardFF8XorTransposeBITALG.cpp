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
    defined(LEO_FF8XOR_BUILD_AVX512_BITALG_TRANSPOSE)

#include <immintrin.h>
#include <string.h>

namespace leopard { namespace ff8xor { namespace transpose { namespace avx512 {


alignas(64) static const uint8_t kBitGatherControl[64] = {
    0, 8, 16, 24, 32, 40, 48, 56,
    0, 8, 16, 24, 32, 40, 48, 56,
    0, 8, 16, 24, 32, 40, 48, 56,
    0, 8, 16, 24, 32, 40, 48, 56,
    0, 8, 16, 24, 32, 40, 48, 56,
    0, 8, 16, 24, 32, 40, 48, 56,
    0, 8, 16, 24, 32, 40, 48, 56,
    0, 8, 16, 24, 32, 40, 48, 56
};

bool BitalgKernelsBuilt()
{
    return true;
}

uint64_t PackedToPlaneBitalg(
    const void* packed_void,
    void* plane_void,
    uint64_t plane_bytes,
    uint64_t group)
{
    const uint8_t* packed = static_cast<const uint8_t*>(packed_void);
    uint8_t* plane = static_cast<uint8_t*>(plane_void);

    if (group > plane_bytes)
        return group;

    const __m512i base = _mm512_load_si512(kBitGatherControl);
    const __m512i one = _mm512_set1_epi8(1);

    while (plane_bytes - group >= 8)
    {
        const __m512i input = _mm512_loadu_si512(packed + group * 8);
        __m512i control = base;

#define LEO_STORE_PLANE(Plane) do { \
        const __mmask64 bits = \
            _mm512_bitshuffle_epi64_mask(input, control); \
        const uint64_t word = static_cast<uint64_t>(bits); \
        memcpy(plane + plane_bytes * (Plane) + group, \
            &word, sizeof(word)); \
        control = _mm512_add_epi8(control, one); \
    } while (0)

        LEO_STORE_PLANE(0);
        LEO_STORE_PLANE(1);
        LEO_STORE_PLANE(2);
        LEO_STORE_PLANE(3);
        LEO_STORE_PLANE(4);
        LEO_STORE_PLANE(5);
        LEO_STORE_PLANE(6);
        LEO_STORE_PLANE(7);

#undef LEO_STORE_PLANE

        group += 8;
    }
    return group;
}


}}}} // namespace leopard::ff8xor::transpose::avx512

#else

namespace leopard { namespace ff8xor { namespace transpose { namespace avx512 {

bool BitalgKernelsBuilt()
{
    return false;
}

uint64_t PackedToPlaneBitalg(
    const void*, void*, uint64_t, uint64_t start_group)
{
    return start_group;
}

}}}} // namespace leopard::ff8xor::transpose::avx512

#endif
