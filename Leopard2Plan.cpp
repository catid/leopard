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

#include "Leopard2Plan.h"

#include <string.h>

namespace leopard2_internal {

uint8_t Log2PowerOfTwo(uint32_t size)
{
    if (size < 2 || (size & (size - 1)) != 0)
        return 0;
    uint8_t result = 0;
    while (size > 1)
    {
        size >>= 1;
        ++result;
    }
    return result;
}

size_t OutputDependencyBitCount(uint32_t transform_size)
{
    const uint8_t log2_size = Log2PowerOfTwo(transform_size);
    if (log2_size == 0)
        return 0;

    size_t result = 0;
    unsigned mip_level = log2_size;
    while (mip_level >= 2)
    {
        result += transform_size >> mip_level;
        mip_level -= 2;
    }
    if (mip_level == 1)
        result += transform_size >> 1;
    return result;
}

size_t OutputDependencyWordCount(uint32_t transform_size)
{
    const size_t bits = OutputDependencyBitCount(transform_size);
    return (bits + 63u) / 64u;
}

bool BuildOutputDependencies(
    uint32_t transform_size,
    const uint32_t* requested_coordinates,
    size_t requested_count,
    uint64_t* words,
    size_t word_count)
{
    const uint8_t log2_size = Log2PowerOfTwo(transform_size);
    const size_t expected_words = OutputDependencyWordCount(transform_size);
    if (log2_size == 0 || word_count != expected_words ||
        (word_count != 0 && !words) ||
        (requested_count != 0 && !requested_coordinates))
        return false;

    // Validate the complete request before publishing any part of the
    // schedule.  Callers may retain and reuse their previous immutable
    // schedule when construction fails, so malformed coordinates must not
    // leave a partially cleared or partially rebuilt bitmap behind.
    for (size_t requested = 0; requested < requested_count; ++requested)
        if (requested_coordinates[requested] >= transform_size)
            return false;

    memset(words, 0, word_count * sizeof(uint64_t));
    for (size_t requested = 0; requested < requested_count; ++requested)
    {
        const uint32_t coordinate = requested_coordinates[requested];
        uint32_t groups_before = 0;
        unsigned mip_level = log2_size;
        while (mip_level >= 2)
        {
            const uint32_t index = groups_before + (coordinate >> mip_level);
            words[index >> 6] |= static_cast<uint64_t>(1) << (index & 63u);
            groups_before += transform_size >> mip_level;
            mip_level -= 2;
        }
        if (mip_level == 1)
        {
            const uint32_t index = groups_before + (coordinate >> 1);
            words[index >> 6] |= static_cast<uint64_t>(1) << (index & 63u);
        }
    }
    return true;
}

} // namespace leopard2_internal
