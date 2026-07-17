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
#include <vector>

namespace leopard { namespace backend { struct Ops; }}

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

enum PrunedTransformOperationFlags
{
    PrunedLiveX = 1u << 0,
    PrunedLiveY = 1u << 1,
    PrunedNeedX = 1u << 2,
    PrunedNeedY = 1u << 3,
    PrunedWriteX = 1u << 4,
    PrunedWriteY = 1u << 5
};

struct PrunedTransformOperation
{
    uint32_t x;
    uint32_t y;
    uint16_t multiplier_log;
    uint8_t flags;
};

/*
    Immutable flat schedule for one exact parent-preserving LCH transform.
    input_mask marks the only coordinates allowed to be nonzero at execution;
    output_mask marks the coordinates whose final values are observable.  The
    executor may modify dead coordinates.  Plan construction is setup work and
    execution performs no allocation or mask branching.
*/
struct PrunedTransformPlan
{
    uint32_t size;
    uint32_t shift;
    uint16_t zero_multiplier_log;
    bool inverse;
    std::vector<uint8_t> input_mask;
    std::vector<uint8_t> output_mask;
    std::vector<PrunedTransformOperation> operations;
    // Requested coordinates that structural analysis proves are zero.  An
    // executor may use dead slots as temporary storage, then clears these at
    // final scatter so requested zero outputs remain exact.
    std::vector<uint32_t> zero_outputs;
    size_t full_butterfly_count;
    size_t one_output_butterflies;
    size_t input_zero_specializations;
    size_t zero_multiplier_butterflies;
    size_t one_multiplier_butterflies;

    PrunedTransformPlan()
        : size(0)
        , shift(0)
        , zero_multiplier_log(0)
        , inverse(false)
        , full_butterfly_count(0)
        , one_output_butterflies(0)
        , input_zero_specializations(0)
        , zero_multiplier_butterflies(0)
        , one_multiplier_butterflies(0)
    {}
};

typedef uint16_t (*PrunedMultiplierLogProvider)(
    const void* context,
    uint32_t storage_index);

// Builds into temporary storage and publishes only on success.  size and
// field_order must be powers of two, shift must name an aligned in-field coset,
// and both masks contain exactly size bytes with values in {0,1}.
bool CompilePrunedTransformPlan(
    uint32_t field_order,
    uint16_t zero_multiplier_log,
    uint32_t size,
    uint32_t shift,
    bool inverse,
    const uint8_t* input_mask,
    const uint8_t* output_mask,
    PrunedMultiplierLogProvider multiplier_log,
    const void* multiplier_context,
    PrunedTransformPlan& plan);

// Executes a trusted immutable schedule without allocating.  The caller must
// provide size pairwise-disjoint shard buffers, keep every masked-off input at
// mathematical zero, and use complete GF16 symbols for a GF16 plan.  Dead
// coordinates may be used as temporary storage and are not preserved.
bool ExecutePrunedTransformPlan(
    const leopard::backend::Ops& ops,
    uint64_t byte_count,
    const PrunedTransformPlan& plan,
    void** work);

uint8_t Log2PowerOfTwo(uint32_t size);
size_t OutputDependencyBitCount(uint32_t transform_size);
size_t OutputDependencyWordCount(uint32_t transform_size);

// transform_size must be a power of two in [2, 2^31].  Returns false for
// malformed sizes, coordinates, pointers, or storage lengths without changing
// caller storage.  In particular, transform_size == 1 is not a valid schedule.
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
