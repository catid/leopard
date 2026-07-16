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

#include <stddef.h>
#include <stdint.h>

namespace leopard2_internal {

/*
    Compact immutable dependency schedule for Leopard's fused-two-layer FFT.

    One bit is stored for every group queried by the traversal: groups of N,
    N/4, N/16, ... coordinates, followed by groups of two when log2(N) is odd.
    Levels are concatenated without per-level padding.  This representation is
    field-neutral: N=256 requires 85 bits (16 bytes), and N=65536 requires
    21845 bits (2736 bytes).
*/
struct OutputDependencyView
{
    const uint64_t* words;
    uint32_t word_count;
    uint8_t log2_size;

    bool IsNeeded(unsigned mip_level, unsigned coordinate) const
    {
        if (!words || log2_size >= 32 ||
            coordinate >= (static_cast<uint32_t>(1) << log2_size) ||
            mip_level == 0 || mip_level > log2_size ||
            ((static_cast<unsigned>(log2_size) - mip_level) & 1u) != 0)
            return false;
        const unsigned level = (static_cast<unsigned>(log2_size) - mip_level) / 2;
        const uint32_t groups_before =
            ((static_cast<uint32_t>(1) << (level * 2)) - 1u) / 3u;
        const uint32_t index = groups_before + (coordinate >> mip_level);
        return index < word_count * 64u &&
            0 != (words[index >> 6] & (static_cast<uint64_t>(1) << (index & 63u)));
    }
};

struct DecodeOutputBlock
{
    uint32_t block;
    uint32_t requested_prefix;
    uint32_t requested_begin;
    uint32_t requested_end;
};

uint8_t Log2PowerOfTwo(uint32_t size);
size_t OutputDependencyBitCount(uint32_t transform_size);
size_t OutputDependencyWordCount(uint32_t transform_size);

// Returns false for malformed sizes, coordinates, or storage lengths.
bool BuildOutputDependencies(
    uint32_t transform_size,
    const uint32_t* requested_coordinates,
    size_t requested_count,
    uint64_t* words,
    size_t word_count);

inline OutputDependencyView MakeOutputDependencyView(
    uint32_t transform_size,
    const uint64_t* words,
    size_t word_count)
{
    OutputDependencyView view = {
        words,
        static_cast<uint32_t>(word_count),
        Log2PowerOfTwo(transform_size)
    };
    return view;
}

} // namespace leopard2_internal
