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
    // Forward plans retain the root-layer multiplier so execution can begin
    // with an immutable-source out-of-place butterfly instead of copying the
    // complete input block to a second workspace.  The inverse accumulating
    // executor uses the same multiplier for its deferred root boundary;
    // materialized inverse execution already carries it in root operations.
    uint16_t first_layer_multiplier_log;
    bool inverse;
    std::vector<uint8_t> input_mask;
    std::vector<uint8_t> output_mask;
    std::vector<PrunedTransformOperation> operations;
    // Liveness/dependency flags for the complete inverse root layer, including
    // identity rows omitted from the compact in-place operation list.  Root
    // coordinates and the multiplier are implicit in the plan identity, so
    // one byte per pair is sufficient.  An accumulating executor consumes
    // this immutable boundary as its final sink instead of storing and
    // rereading a complete syndrome block.  Forward plans leave it empty.
    std::vector<uint8_t> inverse_accumulation_flags;
    // Each sorted index replaces operations i..i+3 with one mature fused
    // four-way backend call.  The compiler emits these descriptors only for a
    // complete two-layer subtransform whose four inputs and four outputs are
    // live; the radix-2 entries remain inspectable without a tagged union.
    std::vector<uint32_t> fused_four_starts;
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
        , first_layer_multiplier_log(0)
        , inverse(false)
        , full_butterfly_count(0)
        , one_output_butterflies(0)
        , input_zero_specializations(0)
        , zero_multiplier_butterflies(0)
        , one_multiplier_butterflies(0)
    {}
};

// One block-local transform schedule in a larger active parent.  Entries are
// sorted by block so Algorithm 5 execution can consume them without lookup
// allocation or branches in the byte inner loop.
struct PrunedTransformBlock
{
    uint32_t block;
    PrunedTransformPlan plan;
};

typedef uint16_t (*PrunedMultiplierLogProvider)(
    const void* context,
    uint32_t storage_index);

/*
    A call-local, allocation-free forward-transform schedule.  A two-bit
    output mask corresponds to each radix-2 butterfly in the ordinary LCH
    traversal, so a retained one-row operation does not write its dead peer.
    The schedule is compiled into caller scratch from the exact requested
    output mask, then reused for every byte tile in the encode call.  It is
    intentionally a non-owning view: persistent decode plans continue to use
    PrunedTransformPlan above, while encode request masks need not mutate an
    immutable codec or allocate a heap object in the hot path.
*/
struct SparseForwardPlanStats
{
    size_t full_butterfly_count;
    size_t retained_butterfly_count;
    size_t one_output_butterflies;
    size_t fused_four_groups;

    SparseForwardPlanStats()
        : full_butterfly_count(0)
        , retained_butterfly_count(0)
        , one_output_butterflies(0)
        , fused_four_groups(0)
    {}
};

struct SparseForwardPlanBatchView
{
    const uint8_t* operation_masks;
    size_t operation_stride;
    uint32_t block_count;
};

size_t SparseForwardButterflyCount(uint32_t transform_size);
size_t SparseForwardRetainedBytes(uint32_t transform_size);
size_t SparseForwardDependencyBytes(uint32_t transform_size);
size_t CountSparseForwardRetainedButterflies(
    uint32_t transform_size,
    const uint8_t* operation_masks,
    size_t retained_bytes);
size_t PrefixForwardButterflyCount(
    uint32_t transform_size,
    uint32_t requested_prefix);

// dependency_workspace initially contains a packed transform_size-bit output
// mask and is consumed in place.  operation_masks is cleared before
// publication and stores two output bits per butterfly.  Failure leaves stats
// empty; callers discard scratch contents.
bool CompileSparseForwardPlan(
    uint32_t field_order,
    uint16_t zero_multiplier_log,
    uint32_t transform_size,
    uint32_t shift,
    uint8_t* dependency_workspace,
    size_t dependency_bytes,
    uint8_t* operation_masks,
    size_t retained_bytes,
    PrunedMultiplierLogProvider multiplier_log,
    const void* multiplier_context,
    SparseForwardPlanStats& stats);

// Executes a trusted call-local schedule without allocation.  The in-place
// form consumes a complete coefficient block.  The source form keeps that
// block immutable and materializes only retained root butterflies in work.
bool ExecuteSparseForwardPlan(
    const leopard::backend::Ops& ops,
    uint64_t byte_count,
    uint32_t transform_size,
    uint32_t shift,
    uint16_t zero_multiplier_log,
    const uint8_t* operation_masks,
    size_t retained_bytes,
    PrunedMultiplierLogProvider multiplier_log,
    const void* multiplier_context,
    void** work);

bool ExecuteSparseForwardPlanFromSources(
    const leopard::backend::Ops& ops,
    uint64_t byte_count,
    uint32_t transform_size,
    uint32_t shift,
    uint16_t zero_multiplier_log,
    const uint8_t* operation_masks,
    size_t retained_bytes,
    PrunedMultiplierLogProvider multiplier_log,
    const void* multiplier_context,
    void* const* source,
    void** work);

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

// Inverse-only companion that executes every non-root operation in work and
// XORs the exact requested root outputs into accumulator.  The source work
// buffers are dead after the call.  work[0..size) and
// accumulator[0..size) must each be pairwise disjoint, and no work byte may
// overlap an accumulator byte.  A rejected plan returns false before changing
// either array; caller-owned materialize-then-XOR therefore remains a safe
// fallback.  Masked-off outputs leave the corresponding accumulator unchanged.
bool ExecutePrunedInverseTransformPlanAccumulate(
    const leopard::backend::Ops& ops,
    uint64_t byte_count,
    const PrunedTransformPlan& plan,
    void** work,
    void** accumulator);

// Forward-only companion for a plan whose input mask is completely live.
// source is immutable and pairwise disjoint from every output buffer in work.
// The root radix-2 layer is executed directly out of place, and the exact
// sparse schedule resumes in work without a whole-size input copy.  This is
// used by high-rate Algorithm 5 to evaluate one reusable coefficient block on
// multiple requested message cosets.
bool ExecutePrunedForwardTransformPlanFromSources(
    const leopard::backend::Ops& ops,
    uint64_t byte_count,
    const PrunedTransformPlan& plan,
    void* const* source,
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
