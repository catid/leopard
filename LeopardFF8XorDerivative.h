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

#pragma once

#include "LeopardCommon.h"

#ifdef LEO_HAS_FF8

namespace leopard { namespace ff8xor { namespace derivative_detail {


// The formal-derivative/top-FFT boundary has only eight source arities.  Each
// specialization keeps every payload value statically named: Extra0..Extra6
// are selected before entering the vector/word loop and there is no indexed
// register array or gate interpreter in the loop.
template <typename Ops>
struct RowChunk
{
    typedef typename Ops::Value Value;

    static LEO_FORCE_INLINE void StoreResults(
        uint8_t* left,
        uint8_t* right,
        uint64_t offset,
        Value a,
        Value r0,
        Value y)
    {
        const Value x = Ops::Xor(a, r0);
        Ops::Store(left + offset, x);
        Ops::Store(right + offset, y);
    }

    static LEO_FORCE_INLINE void Apply0(
        uint8_t* left,
        uint8_t* right,
        uint64_t offset)
    {
        const Value a = Ops::Load(left + offset);
        const Value r0 = Ops::Load(right + offset);
        StoreResults(left, right, offset, a, r0, a);
    }

    static LEO_FORCE_INLINE void Apply1(
        uint8_t* left,
        uint8_t* right,
        const uint8_t* extra0,
        uint64_t offset)
    {
        const Value a = Ops::Load(left + offset);
        const Value r0 = Ops::Load(right + offset);
        const Value s0 = Ops::Load(extra0 + offset);
        const Value y = Ops::Xor(a, s0);
        StoreResults(left, right, offset, a, r0, y);
    }

    static LEO_FORCE_INLINE void Apply2(
        uint8_t* left,
        uint8_t* right,
        const uint8_t* extra0,
        const uint8_t* extra1,
        uint64_t offset)
    {
        const Value a = Ops::Load(left + offset);
        const Value r0 = Ops::Load(right + offset);
        const Value s0 = Ops::Load(extra0 + offset);
        const Value s1 = Ops::Load(extra1 + offset);
        const Value y = Ops::Xor3(a, s0, s1);
        StoreResults(left, right, offset, a, r0, y);
    }

    static LEO_FORCE_INLINE void Apply3(
        uint8_t* left,
        uint8_t* right,
        const uint8_t* extra0,
        const uint8_t* extra1,
        const uint8_t* extra2,
        uint64_t offset)
    {
        const Value a = Ops::Load(left + offset);
        const Value r0 = Ops::Load(right + offset);
        const Value s0 = Ops::Load(extra0 + offset);
        const Value s1 = Ops::Load(extra1 + offset);
        const Value s2 = Ops::Load(extra2 + offset);
        Value y = Ops::Xor3(a, s0, s1);
        y = Ops::Xor(y, s2);
        StoreResults(left, right, offset, a, r0, y);
    }

    static LEO_FORCE_INLINE void Apply4(
        uint8_t* left,
        uint8_t* right,
        const uint8_t* extra0,
        const uint8_t* extra1,
        const uint8_t* extra2,
        const uint8_t* extra3,
        uint64_t offset)
    {
        const Value a = Ops::Load(left + offset);
        const Value r0 = Ops::Load(right + offset);
        const Value s0 = Ops::Load(extra0 + offset);
        const Value s1 = Ops::Load(extra1 + offset);
        const Value s2 = Ops::Load(extra2 + offset);
        const Value s3 = Ops::Load(extra3 + offset);
        Value y = Ops::Xor3(a, s0, s1);
        y = Ops::Xor3(y, s2, s3);
        StoreResults(left, right, offset, a, r0, y);
    }

    static LEO_FORCE_INLINE void Apply5(
        uint8_t* left,
        uint8_t* right,
        const uint8_t* extra0,
        const uint8_t* extra1,
        const uint8_t* extra2,
        const uint8_t* extra3,
        const uint8_t* extra4,
        uint64_t offset)
    {
        const Value a = Ops::Load(left + offset);
        const Value r0 = Ops::Load(right + offset);
        const Value s0 = Ops::Load(extra0 + offset);
        const Value s1 = Ops::Load(extra1 + offset);
        const Value s2 = Ops::Load(extra2 + offset);
        const Value s3 = Ops::Load(extra3 + offset);
        const Value s4 = Ops::Load(extra4 + offset);
        Value y = Ops::Xor3(a, s0, s1);
        y = Ops::Xor3(y, s2, s3);
        y = Ops::Xor(y, s4);
        StoreResults(left, right, offset, a, r0, y);
    }

    static LEO_FORCE_INLINE void Apply6(
        uint8_t* left,
        uint8_t* right,
        const uint8_t* extra0,
        const uint8_t* extra1,
        const uint8_t* extra2,
        const uint8_t* extra3,
        const uint8_t* extra4,
        const uint8_t* extra5,
        uint64_t offset)
    {
        const Value a = Ops::Load(left + offset);
        const Value r0 = Ops::Load(right + offset);
        const Value s0 = Ops::Load(extra0 + offset);
        const Value s1 = Ops::Load(extra1 + offset);
        const Value s2 = Ops::Load(extra2 + offset);
        const Value s3 = Ops::Load(extra3 + offset);
        const Value s4 = Ops::Load(extra4 + offset);
        const Value s5 = Ops::Load(extra5 + offset);
        Value y = Ops::Xor3(a, s0, s1);
        y = Ops::Xor3(y, s2, s3);
        y = Ops::Xor3(y, s4, s5);
        StoreResults(left, right, offset, a, r0, y);
    }

    static LEO_FORCE_INLINE void Apply7(
        uint8_t* left,
        uint8_t* right,
        const uint8_t* extra0,
        const uint8_t* extra1,
        const uint8_t* extra2,
        const uint8_t* extra3,
        const uint8_t* extra4,
        const uint8_t* extra5,
        const uint8_t* extra6,
        uint64_t offset)
    {
        const Value a = Ops::Load(left + offset);
        const Value r0 = Ops::Load(right + offset);
        const Value s0 = Ops::Load(extra0 + offset);
        const Value s1 = Ops::Load(extra1 + offset);
        const Value s2 = Ops::Load(extra2 + offset);
        const Value s3 = Ops::Load(extra3 + offset);
        const Value s4 = Ops::Load(extra4 + offset);
        const Value s5 = Ops::Load(extra5 + offset);
        const Value s6 = Ops::Load(extra6 + offset);
        Value y = Ops::Xor3(a, s0, s1);
        y = Ops::Xor3(y, s2, s3);
        y = Ops::Xor3(y, s4, s5);
        y = Ops::Xor(y, s6);
        StoreResults(left, right, offset, a, r0, y);
    }
};

template <typename Ops>
uint64_t ApplyRow(
    unsigned extra_count,
    void* left_void,
    void* right_void,
    const void* extra0_void,
    const void* extra1_void,
    const void* extra2_void,
    const void* extra3_void,
    const void* extra4_void,
    const void* extra5_void,
    const void* extra6_void,
    uint64_t buffer_bytes,
    uint64_t offset)
{
    uint8_t* left = reinterpret_cast<uint8_t*>(left_void);
    uint8_t* right = reinterpret_cast<uint8_t*>(right_void);
    const uint8_t* extra0 = reinterpret_cast<const uint8_t*>(extra0_void);
    const uint8_t* extra1 = reinterpret_cast<const uint8_t*>(extra1_void);
    const uint8_t* extra2 = reinterpret_cast<const uint8_t*>(extra2_void);
    const uint8_t* extra3 = reinterpret_cast<const uint8_t*>(extra3_void);
    const uint8_t* extra4 = reinterpret_cast<const uint8_t*>(extra4_void);
    const uint8_t* extra5 = reinterpret_cast<const uint8_t*>(extra5_void);
    const uint8_t* extra6 = reinterpret_cast<const uint8_t*>(extra6_void);

    switch (extra_count)
    {
    case 0:
        while (buffer_bytes - offset >= Ops::kBytes)
        {
            RowChunk<Ops>::Apply0(left, right, offset);
            offset += Ops::kBytes;
        }
        break;
    case 1:
        while (buffer_bytes - offset >= Ops::kBytes)
        {
            RowChunk<Ops>::Apply1(left, right, extra0, offset);
            offset += Ops::kBytes;
        }
        break;
    case 2:
        while (buffer_bytes - offset >= Ops::kBytes)
        {
            RowChunk<Ops>::Apply2(
                left, right, extra0, extra1, offset);
            offset += Ops::kBytes;
        }
        break;
    case 3:
        while (buffer_bytes - offset >= Ops::kBytes)
        {
            RowChunk<Ops>::Apply3(
                left, right, extra0, extra1, extra2, offset);
            offset += Ops::kBytes;
        }
        break;
    case 4:
        while (buffer_bytes - offset >= Ops::kBytes)
        {
            RowChunk<Ops>::Apply4(
                left, right, extra0, extra1, extra2, extra3, offset);
            offset += Ops::kBytes;
        }
        break;
    case 5:
        while (buffer_bytes - offset >= Ops::kBytes)
        {
            RowChunk<Ops>::Apply5(
                left, right, extra0, extra1, extra2, extra3, extra4,
                offset);
            offset += Ops::kBytes;
        }
        break;
    case 6:
        while (buffer_bytes - offset >= Ops::kBytes)
        {
            RowChunk<Ops>::Apply6(
                left, right, extra0, extra1, extra2, extra3, extra4,
                extra5, offset);
            offset += Ops::kBytes;
        }
        break;
    case 7:
        while (buffer_bytes - offset >= Ops::kBytes)
        {
            RowChunk<Ops>::Apply7(
                left, right, extra0, extra1, extra2, extra3, extra4,
                extra5, extra6, offset);
            offset += Ops::kBytes;
        }
        break;
    default:
        break;
    }

    return offset;
}


}}} // namespace leopard::ff8xor::derivative_detail

#endif // LEO_HAS_FF8
