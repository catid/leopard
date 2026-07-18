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
#include "Leopard2Backend.h"

#include <string.h>
#include <utility>

namespace {

struct RawPrunedOperation
{
    uint32_t x;
    uint32_t y;
    uint16_t multiplier_log;
};

static bool IsPowerOfTwo(uint32_t value)
{
    return value != 0 && (value & (value - 1)) == 0;
}

static void BuildRawPrunedOperations(
    uint32_t start,
    uint32_t size,
    uint32_t coset,
    bool inverse,
    leopard2_internal::PrunedMultiplierLogProvider multiplier_log,
    const void* multiplier_context,
    std::vector<RawPrunedOperation>& operations)
{
    if (size == 1)
        return;

    const uint32_t half = size >> 1;
    if (size == 2)
    {
        const uint16_t log_m = multiplier_log(
            multiplier_context, coset + half);
        const RawPrunedOperation operation = {
            start,
            start + 1,
            log_m
        };
        operations.push_back(operation);
        return;
    }

    // Group every pair of transform layers in the exact order consumed by
    // Butterfly4.  Operations for different offsets and disjoint child
    // subspaces commute, so this is the same radix-2 DAG as the historical
    // root/child traversal with complete two-layer regions made contiguous.
    const uint32_t quarter = size >> 2;
    const uint16_t log_m01 = multiplier_log(
        multiplier_context, coset + quarter);
    const uint16_t log_m23 = multiplier_log(
        multiplier_context, coset + half + quarter);
    const uint16_t log_m02 = multiplier_log(
        multiplier_context, coset + half);

    if (inverse)
    {
        BuildRawPrunedOperations(
            start, quarter, coset, inverse,
            multiplier_log, multiplier_context, operations);
        BuildRawPrunedOperations(
            start + quarter, quarter, coset + quarter, inverse,
            multiplier_log, multiplier_context, operations);
        BuildRawPrunedOperations(
            start + half, quarter, coset + half, inverse,
            multiplier_log, multiplier_context, operations);
        BuildRawPrunedOperations(
            start + half + quarter, quarter, coset + half + quarter,
            inverse, multiplier_log, multiplier_context, operations);
    }

    for (uint32_t offset = 0; offset < quarter; ++offset)
    {
        const uint32_t value0 = start + offset;
        const uint32_t value1 = value0 + quarter;
        const uint32_t value2 = value0 + half;
        const uint32_t value3 = value2 + quarter;
        if (inverse)
        {
            const RawPrunedOperation grouped[] = {
                { value0, value1, log_m01 },
                { value2, value3, log_m23 },
                { value0, value2, log_m02 },
                { value1, value3, log_m02 }
            };
            operations.insert(operations.end(), grouped, grouped + 4);
        }
        else
        {
            const RawPrunedOperation grouped[] = {
                { value0, value2, log_m02 },
                { value1, value3, log_m02 },
                { value0, value1, log_m01 },
                { value2, value3, log_m23 }
            };
            operations.insert(operations.end(), grouped, grouped + 4);
        }
    }

    if (!inverse)
    {
        BuildRawPrunedOperations(
            start, quarter, coset, inverse,
            multiplier_log, multiplier_context, operations);
        BuildRawPrunedOperations(
            start + quarter, quarter, coset + quarter, inverse,
            multiplier_log, multiplier_context, operations);
        BuildRawPrunedOperations(
            start + half, quarter, coset + half, inverse,
            multiplier_log, multiplier_context, operations);
        BuildRawPrunedOperations(
            start + half + quarter, quarter, coset + half + quarter,
            inverse, multiplier_log, multiplier_context, operations);
    }
}

static bool CoefficientAIsNonzero(bool inverse, bool multiplier_one)
{
    return !inverse || !multiplier_one;
}

static bool CoefficientAIsOne(bool inverse, bool multiplier_zero)
{
    return !inverse || multiplier_zero;
}

static bool CoefficientDIsNonzero(bool inverse, bool multiplier_one)
{
    return inverse || !multiplier_one;
}

static bool CoefficientDIsOne(bool inverse, bool multiplier_zero)
{
    return inverse || multiplier_zero;
}

struct PrunedExecutionOps
{
    leopard::backend::FixedMultiply multiply;
    leopard::backend::FixedMultiply multiply_add;
    leopard::backend::Butterfly2 butterfly;
    leopard::backend::Butterfly4 butterfly_four;
    leopard::backend::FFTButterfly2Out butterfly_out;
};

struct FusedFourMatch
{
    uint32_t value0;
    uint32_t value1;
    uint32_t value2;
    uint32_t value3;
    uint16_t multiplier_log01;
    uint16_t multiplier_log23;
    uint16_t multiplier_log02;
};

static bool MatchFusedFour(
    const std::vector<leopard2_internal::PrunedTransformOperation>& operations,
    size_t start,
    bool inverse,
    FusedFourMatch& match)
{
    if (start > operations.size() || operations.size() - start < 4)
        return false;

    const leopard2_internal::PrunedTransformOperation& first =
        operations[start];
    const leopard2_internal::PrunedTransformOperation& second =
        operations[start + 1];
    const leopard2_internal::PrunedTransformOperation& third =
        operations[start + 2];
    const leopard2_internal::PrunedTransformOperation& fourth =
        operations[start + 3];
    const uint8_t complete =
        leopard2_internal::PrunedLiveX |
        leopard2_internal::PrunedLiveY |
        leopard2_internal::PrunedNeedX |
        leopard2_internal::PrunedNeedY;
    if ((first.flags & complete) != complete ||
        (second.flags & complete) != complete ||
        (third.flags & complete) != complete ||
        (fourth.flags & complete) != complete)
        return false;

    if (inverse)
    {
        match.value0 = first.x;
        match.value1 = first.y;
        match.value2 = second.x;
        match.value3 = second.y;
        match.multiplier_log01 = first.multiplier_log;
        match.multiplier_log23 = second.multiplier_log;
        match.multiplier_log02 = third.multiplier_log;
        if (third.x != match.value0 || third.y != match.value2 ||
            fourth.x != match.value1 || fourth.y != match.value3 ||
            fourth.multiplier_log != match.multiplier_log02)
            return false;
    }
    else
    {
        match.value0 = first.x;
        match.value2 = first.y;
        match.value1 = second.x;
        match.value3 = second.y;
        match.multiplier_log02 = first.multiplier_log;
        match.multiplier_log01 = third.multiplier_log;
        match.multiplier_log23 = fourth.multiplier_log;
        if (second.multiplier_log != match.multiplier_log02 ||
            third.x != match.value0 || third.y != match.value1 ||
            fourth.x != match.value2 || fourth.y != match.value3)
            return false;
    }

    return match.value0 != match.value1 &&
        match.value0 != match.value2 &&
        match.value0 != match.value3 &&
        match.value1 != match.value2 &&
        match.value1 != match.value3 &&
        match.value2 != match.value3;
}

static bool ExecutePrunedOperation(
    const leopard::backend::Ops& ops,
    const PrunedExecutionOps& selected,
    bool inverse,
    uint16_t zero_multiplier_log,
    uint64_t byte_count,
    const leopard2_internal::PrunedTransformOperation& operation,
    void* x,
    void* y)
{
    const uint8_t flags = operation.flags;
    const bool live_x = 0 != (flags & leopard2_internal::PrunedLiveX);
    const bool live_y = 0 != (flags & leopard2_internal::PrunedLiveY);
    const bool write_x = 0 != (flags & leopard2_internal::PrunedWriteX);
    const bool write_y = 0 != (flags & leopard2_internal::PrunedWriteY);
    const uint16_t log_m = operation.multiplier_log;
    const bool multiplier_zero = log_m == zero_multiplier_log;
    const bool multiplier_one = log_m == 0;

    if (live_x && live_y)
    {
        // A retained row may leave its peer dead, but running the complete
        // two-way butterfly is still exact and lets the mature backend reuse
        // its one fixed product.  Zero and one avoid table lookup entirely.
        if (multiplier_zero)
        {
            ops.xor_memory(y, x, byte_count);
        }
        else if (multiplier_one)
        {
            if (inverse)
            {
                ops.xor_memory(y, x, byte_count);
                ops.xor_memory(x, y, byte_count);
            }
            else
            {
                ops.xor_memory(x, y, byte_count);
                ops.xor_memory(y, x, byte_count);
            }
        }
        else
        {
            selected.butterfly(x, y, log_m, byte_count);
        }
        return true;
    }

    if (live_x)
    {
        if (!inverse)
        {
            // [x, 0] -> [x, x].
            if (write_y)
                memcpy(y, x, static_cast<size_t>(byte_count));
            return true;
        }

        // [x, 0] -> [(m + 1) x, x].  The dead y slot safely preserves x
        // while the in-place multiply-add updates the first row.
        if (write_y || (write_x && !multiplier_zero))
            memcpy(y, x, static_cast<size_t>(byte_count));
        if (write_x && !multiplier_zero)
            selected.multiply_add(x, y, log_m, byte_count);
        return true;
    }

    if (live_y)
    {
        // Forward: [0, y] -> [m y, (m + 1) y].
        // Inverse: [0, y] -> [m y, y].
        // If only the second forward row is needed, x remains available as a
        // dead temporary so the direct matrix row still needs one product.
        if (write_x || (!inverse && write_y))
        {
            if (multiplier_zero)
                memset(x, 0, static_cast<size_t>(byte_count));
            else if (multiplier_one)
                memcpy(x, y, static_cast<size_t>(byte_count));
            else
                selected.multiply(x, y, log_m, byte_count);
        }
        if (!inverse && write_y)
            ops.xor_memory(y, x, byte_count);
        return true;
    }

    return !write_x && !write_y;
}

static size_t SparseButterflyCountUnchecked(uint32_t size)
{
    size_t log2_size = 0;
    for (uint32_t value = size; value > 1; value >>= 1)
        ++log2_size;
    return static_cast<size_t>(size >> 1) * log2_size;
}

static bool PackedBit(
    const uint8_t* bits,
    size_t bit_count,
    size_t index)
{
    return index < bit_count &&
        0 != (bits[index >> 3] & static_cast<uint8_t>(1u << (index & 7u)));
}

static void AssignPackedBit(uint8_t* bits, size_t index, bool value)
{
    const uint8_t mask = static_cast<uint8_t>(1u << (index & 7u));
    if (value)
        bits[index >> 3] |= mask;
    else
        bits[index >> 3] &= static_cast<uint8_t>(~mask);
}

static uint8_t SparseOperationMask(
    const uint8_t* masks,
    size_t operation_count,
    size_t index)
{
    if (index >= operation_count)
        return 0;
    return static_cast<uint8_t>(
        (masks[index >> 2] >> ((index & 3U) * 2U)) & 3U);
}

static void SetSparseOperationMask(
    uint8_t* masks,
    size_t index,
    uint8_t value)
{
    const unsigned shift = static_cast<unsigned>((index & 3U) * 2U);
    const uint8_t clear = static_cast<uint8_t>(~(3U << shift));
    masks[index >> 2] = static_cast<uint8_t>(
        (masks[index >> 2] & clear) | (value << shift));
}

struct SparseCompileContext
{
    uint16_t zero_multiplier_log;
    uint8_t* needed;
    size_t needed_count;
    uint8_t* operation_masks;
    size_t operation_index;
    leopard2_internal::PrunedMultiplierLogProvider multiplier_log;
    const void* multiplier_context;
    leopard2_internal::SparseForwardPlanStats* stats;
    bool valid;
};

static uint16_t SparseMultiplier(
    SparseCompileContext& context,
    uint32_t storage_index)
{
    const uint16_t result = context.multiplier_log(
        context.multiplier_context, storage_index);
    if (result > context.zero_multiplier_log)
        context.valid = false;
    return result;
}

static void CompileSparseOperationReverse(
    SparseCompileContext& context,
    uint32_t x,
    uint32_t y,
    uint16_t multiplier_log)
{
    if (!context.valid || context.operation_index == 0)
    {
        context.valid = false;
        return;
    }
    const size_t operation = --context.operation_index;
    const bool need_x = PackedBit(context.needed, context.needed_count, x);
    const bool need_y = PackedBit(context.needed, context.needed_count, y);
    if (!need_x && !need_y)
        return;

    SetSparseOperationMask(context.operation_masks, operation,
        static_cast<uint8_t>((need_x ? 1U : 0U) | (need_y ? 2U : 0U)));
    ++context.stats->retained_butterfly_count;
    if (need_x != need_y)
        ++context.stats->one_output_butterflies;

    const bool multiplier_zero =
        multiplier_log == context.zero_multiplier_log;
    const bool multiplier_one = multiplier_log == 0;
    // Forward butterfly rows are x + m*y and x + (m+1)*y.
    AssignPackedBit(context.needed, x, need_x || need_y);
    AssignPackedBit(context.needed, y,
        (need_x && !multiplier_zero) ||
        (need_y && !multiplier_one));
}

static void CompileSparseNodeReverse(
    SparseCompileContext& context,
    uint32_t start,
    uint32_t size,
    uint32_t coset)
{
    if (!context.valid)
        return;
    if (size == 1)
        return;
    if (size == 2)
    {
        CompileSparseOperationReverse(context, start, start + 1,
            SparseMultiplier(context, coset + 1));
        return;
    }

    const uint32_t quarter = size >> 2;
    const uint32_t half = size >> 1;
    const uint16_t multiplier01 = SparseMultiplier(
        context, coset + quarter);
    const uint16_t multiplier23 = SparseMultiplier(
        context, coset + half + quarter);
    const uint16_t multiplier02 = SparseMultiplier(
        context, coset + half);

    // Reverse the exact forward traversal: four child transforms followed by
    // the grouped two-layer butterflies at this node.
    CompileSparseNodeReverse(context,
        start + half + quarter, quarter, coset + half + quarter);
    CompileSparseNodeReverse(context,
        start + half, quarter, coset + half);
    CompileSparseNodeReverse(context,
        start + quarter, quarter, coset + quarter);
    CompileSparseNodeReverse(context,
        start, quarter, coset);

    for (uint32_t offset = quarter; offset-- > 0;)
    {
        const uint32_t value0 = start + offset;
        const uint32_t value1 = value0 + quarter;
        const uint32_t value2 = value0 + half;
        const uint32_t value3 = value2 + quarter;
        CompileSparseOperationReverse(
            context, value2, value3, multiplier23);
        CompileSparseOperationReverse(
            context, value0, value1, multiplier01);
        CompileSparseOperationReverse(
            context, value1, value3, multiplier02);
        CompileSparseOperationReverse(
            context, value0, value2, multiplier02);
    }
}

struct SparseCountContext
{
    const uint8_t* operation_masks;
    size_t operation_count;
    size_t operation_index;
    size_t fused_four_groups;
};

static void CountSparseFusedNode(
    SparseCountContext& context,
    uint32_t size)
{
    if (size == 1)
        return;
    if (size == 2)
    {
        ++context.operation_index;
        return;
    }
    const uint32_t quarter = size >> 2;
    for (uint32_t offset = 0; offset < quarter; ++offset)
    {
        bool complete = true;
        for (unsigned operation = 0; operation < 4; ++operation)
            complete = complete && SparseOperationMask(
                context.operation_masks, context.operation_count,
                context.operation_index + operation) == 3;
        if (complete)
            ++context.fused_four_groups;
        context.operation_index += 4;
    }
    CountSparseFusedNode(context, quarter);
    CountSparseFusedNode(context, quarter);
    CountSparseFusedNode(context, quarter);
    CountSparseFusedNode(context, quarter);
}

struct SparseExecuteContext
{
    const leopard::backend::Ops* ops;
    PrunedExecutionOps selected;
    uint16_t zero_multiplier_log;
    uint64_t byte_count;
    const uint8_t* operation_masks;
    size_t operation_count;
    size_t operation_index;
    uint32_t skipped_distance;
    leopard2_internal::PrunedMultiplierLogProvider multiplier_log;
    const void* multiplier_context;
    void** work;
    bool valid;
};

static uint16_t SparseMultiplier(
    SparseExecuteContext& context,
    uint32_t storage_index)
{
    const uint16_t result = context.multiplier_log(
        context.multiplier_context, storage_index);
    if (result > context.zero_multiplier_log)
        context.valid = false;
    return result;
}

static bool ExecuteSparseOperation(
    SparseExecuteContext& context,
    size_t operation,
    uint32_t x,
    uint32_t y,
    uint16_t multiplier_log)
{
    const uint8_t output_mask = SparseOperationMask(
        context.operation_masks, context.operation_count, operation);
    if (output_mask == 0 || y - x == context.skipped_distance)
        return true;
    const bool multiplier_zero =
        multiplier_log == context.zero_multiplier_log;
    const bool multiplier_one = multiplier_log == 0;
    void* const x_output = context.work[x];
    void* const y_output = context.work[y];
    if (output_mask == 3)
    {
        if (multiplier_zero)
            context.ops->xor_memory(y_output, x_output, context.byte_count);
        else if (multiplier_one)
        {
            context.ops->xor_memory(x_output, y_output, context.byte_count);
            context.ops->xor_memory(y_output, x_output, context.byte_count);
        }
        else
            context.selected.butterfly(
                x_output, y_output, multiplier_log, context.byte_count);
        return true;
    }
    if (output_mask == 1)
    {
        // x' = x + m*y.  In particular m=0 is an identity and performs no
        // write to either row.
        if (multiplier_zero)
            return true;
        if (multiplier_one)
            context.ops->xor_memory(x_output, y_output, context.byte_count);
        else
            context.selected.multiply_add(
                x_output, y_output, multiplier_log, context.byte_count);
        return true;
    }

    // y' = x + (m+1)*y.  multiply_add supports destination==source in every
    // qualified backend: each vector/scalar lane is loaded before it is stored.
    if (multiplier_zero)
        context.ops->xor_memory(y_output, x_output, context.byte_count);
    else if (multiplier_one)
        memcpy(y_output, x_output, static_cast<size_t>(context.byte_count));
    else
    {
        context.selected.multiply_add(
            y_output, y_output, multiplier_log, context.byte_count);
        context.ops->xor_memory(y_output, x_output, context.byte_count);
    }
    return true;
}

static bool ExecuteSparseOperationFromSources(
    SparseExecuteContext& context,
    size_t operation,
    uint32_t x,
    uint32_t y,
    uint16_t multiplier_log,
    void* const* source)
{
    const uint8_t output_mask = SparseOperationMask(
        context.operation_masks, context.operation_count, operation);
    if (output_mask == 0)
        return true;
    const bool multiplier_zero =
        multiplier_log == context.zero_multiplier_log;
    const bool multiplier_one = multiplier_log == 0;
    const void* const x_source = source[x];
    const void* const y_source = source[y];
    void* const x_output = context.work[x];
    void* const y_output = context.work[y];
    const size_t bytes = static_cast<size_t>(context.byte_count);
    if (output_mask == 3)
    {
        if (!multiplier_zero)
        {
            context.selected.butterfly_out(
                x_source, y_source, x_output, y_output,
                multiplier_log, context.byte_count);
            return true;
        }
        memcpy(x_output, x_source, bytes);
        memcpy(y_output, y_source, bytes);
        context.ops->xor_memory(y_output, x_source, context.byte_count);
        return true;
    }
    if (output_mask == 1)
    {
        if (multiplier_zero)
            memcpy(x_output, x_source, bytes);
        else if (multiplier_one)
        {
            memcpy(x_output, x_source, bytes);
            context.ops->xor_memory(x_output, y_source, context.byte_count);
        }
        else
        {
            context.selected.multiply(
                x_output, y_source, multiplier_log, context.byte_count);
            context.ops->xor_memory(x_output, x_source, context.byte_count);
        }
        return true;
    }

    if (multiplier_one)
        memcpy(y_output, x_source, bytes);
    else
    {
        memcpy(y_output, y_source, bytes);
        if (!multiplier_zero)
            context.selected.multiply_add(
                y_output, y_source, multiplier_log, context.byte_count);
        context.ops->xor_memory(y_output, x_source, context.byte_count);
    }
    return true;
}

static void ExecuteSparseNode(
    SparseExecuteContext& context,
    uint32_t start,
    uint32_t size,
    uint32_t coset)
{
    if (!context.valid)
        return;
    if (size == 1)
        return;
    if (size == 2)
    {
        const size_t operation = context.operation_index++;
        const uint16_t multiplier = SparseMultiplier(context, coset + 1);
        if (context.valid)
            ExecuteSparseOperation(
                context, operation, start, start + 1, multiplier);
        return;
    }

    const uint32_t quarter = size >> 2;
    const uint32_t half = size >> 1;
    const uint16_t multiplier01 = SparseMultiplier(
        context, coset + quarter);
    const uint16_t multiplier23 = SparseMultiplier(
        context, coset + half + quarter);
    const uint16_t multiplier02 = SparseMultiplier(
        context, coset + half);
    if (!context.valid)
        return;
    for (uint32_t offset = 0; offset < quarter; ++offset)
    {
        const uint32_t value0 = start + offset;
        const uint32_t value1 = value0 + quarter;
        const uint32_t value2 = value0 + half;
        const uint32_t value3 = value2 + quarter;
        const size_t first = context.operation_index;
        const bool skipped = half == context.skipped_distance;
        bool complete = !skipped;
        for (unsigned operation = 0; operation < 4; ++operation)
            complete = complete && SparseOperationMask(
                context.operation_masks, context.operation_count,
                first + operation) == 3;
        if (complete)
        {
            context.selected.butterfly_four(
                context.work[value0], context.work[value1],
                context.work[value2], context.work[value3],
                multiplier01, multiplier23, multiplier02,
                context.byte_count);
        }
        else
        {
            ExecuteSparseOperation(context, first, value0, value2,
                multiplier02);
            ExecuteSparseOperation(context, first + 1, value1, value3,
                multiplier02);
            ExecuteSparseOperation(context, first + 2, value0, value1,
                multiplier01);
            ExecuteSparseOperation(context, first + 3, value2, value3,
                multiplier23);
        }
        context.operation_index += 4;
    }

    ExecuteSparseNode(context, start, quarter, coset);
    ExecuteSparseNode(context,
        start + quarter, quarter, coset + quarter);
    ExecuteSparseNode(context,
        start + half, quarter, coset + half);
    ExecuteSparseNode(context,
        start + half + quarter, quarter, coset + half + quarter);
}

static bool ExecuteSparseRootFromSources(
    SparseExecuteContext& context,
    uint32_t size,
    uint32_t shift,
    void* const* source)
{
    const uint32_t half = size >> 1;
    if (size == 2)
    {
        const uint16_t multiplier = SparseMultiplier(context, shift + 1);
        if (!context.valid)
            return false;
        return ExecuteSparseOperationFromSources(
            context, 0, 0, 1, multiplier, source);
    }

    const uint32_t quarter = size >> 2;
    const uint16_t multiplier = SparseMultiplier(context, shift + half);
    if (!context.valid)
        return false;
    for (uint32_t offset = 0; offset < quarter; ++offset)
    {
        const uint32_t pairs[][2] = {
            { offset, offset + half },
            { offset + quarter, offset + half + quarter }
        };
        for (unsigned pair = 0; pair < 2; ++pair)
        {
            const size_t operation = static_cast<size_t>(offset) * 4 + pair;
            const uint32_t x = pairs[pair][0];
            const uint32_t y = pairs[pair][1];
            if (!ExecuteSparseOperationFromSources(
                    context, operation, x, y, multiplier, source))
                return false;
        }
    }
    return context.valid;
}

} // namespace

namespace leopard2_internal {

size_t SparseForwardButterflyCount(uint32_t transform_size)
{
    if (!IsPowerOfTwo(transform_size) || transform_size < 2 ||
        transform_size > 65536U)
        return 0;
    return SparseButterflyCountUnchecked(transform_size);
}

size_t SparseForwardRetainedBytes(uint32_t transform_size)
{
    const size_t butterflies = SparseForwardButterflyCount(transform_size);
    return butterflies == 0 ? 0 : (butterflies + 3U) / 4U;
}

size_t SparseForwardDependencyBytes(uint32_t transform_size)
{
    if (!IsPowerOfTwo(transform_size) || transform_size < 2 ||
        transform_size > 65536U)
        return 0;
    return (static_cast<size_t>(transform_size) + 7U) / 8U;
}

size_t CountSparseForwardRetainedButterflies(
    uint32_t transform_size,
    const uint8_t* operation_masks,
    size_t retained_bytes)
{
    const size_t butterflies = SparseForwardButterflyCount(transform_size);
    if (butterflies == 0 || !operation_masks ||
        retained_bytes != SparseForwardRetainedBytes(transform_size))
        return 0;
    size_t result = 0;
    for (size_t i = 0; i < butterflies; ++i)
        if (SparseOperationMask(operation_masks, butterflies, i) != 0)
            ++result;
    return result;
}

size_t PrefixForwardButterflyCount(
    uint32_t transform_size,
    uint32_t requested_prefix)
{
    if (!IsPowerOfTwo(transform_size) || transform_size < 2 ||
        transform_size > 65536U || requested_prefix > transform_size)
        return 0;
    size_t result = 0;
    uint32_t dist4 = transform_size;
    uint32_t dist = transform_size >> 2;
    for (; dist != 0; dist4 = dist, dist >>= 2)
    {
        for (uint32_t start = 0;
             start < requested_prefix;
             start += dist4)
            result += static_cast<size_t>(dist) * 4U;
    }
    if (dist4 == 2)
        result += (static_cast<size_t>(requested_prefix) + 1U) / 2U;
    return result;
}

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
    SparseForwardPlanStats& stats)
{
    stats = SparseForwardPlanStats();
    const size_t full_count = SparseForwardButterflyCount(transform_size);
    const size_t expected_retained_bytes =
        SparseForwardRetainedBytes(transform_size);
    const size_t expected_dependency_bytes =
        SparseForwardDependencyBytes(transform_size);
    if (!IsPowerOfTwo(field_order) || field_order > 65536U ||
        static_cast<uint32_t>(zero_multiplier_log) + 1U != field_order ||
        full_count == 0 || transform_size > field_order ||
        shift > field_order - transform_size ||
        (shift & (transform_size - 1U)) != 0 ||
        !dependency_workspace || dependency_bytes != expected_dependency_bytes ||
        !operation_masks || retained_bytes != expected_retained_bytes ||
        !multiplier_log || !multiplier_context)
        return false;

    memset(operation_masks, 0, retained_bytes);
    SparseForwardPlanStats candidate;
    candidate.full_butterfly_count = full_count;
    SparseCompileContext context = {
        zero_multiplier_log,
        dependency_workspace,
        transform_size,
        operation_masks,
        full_count,
        multiplier_log,
        multiplier_context,
        &candidate,
        true
    };
    CompileSparseNodeReverse(context, 0, transform_size, shift);
    if (!context.valid || context.operation_index != 0)
        return false;

    SparseCountContext count = {
        operation_masks,
        full_count,
        0,
        0
    };
    CountSparseFusedNode(count, transform_size);
    if (count.operation_index != full_count)
        return false;
    candidate.fused_four_groups = count.fused_four_groups;
    stats = candidate;
    return true;
}

static bool ValidateSparseExecution(
    uint64_t byte_count,
    uint32_t transform_size,
    uint32_t shift,
    uint16_t zero_multiplier_log,
    const uint8_t* operation_masks,
    size_t retained_bytes,
    PrunedMultiplierLogProvider multiplier_log,
    const void* multiplier_context,
    void** work)
{
    const bool gf16 = zero_multiplier_log == 65535U;
    const uint32_t field_order =
        static_cast<uint32_t>(zero_multiplier_log) + 1U;
    if ((zero_multiplier_log != 255U && !gf16) ||
        SparseForwardButterflyCount(transform_size) == 0 ||
        transform_size > field_order ||
        shift > field_order - transform_size ||
        (shift & (transform_size - 1U)) != 0 ||
        retained_bytes != SparseForwardRetainedBytes(transform_size) ||
        !operation_masks || !multiplier_log || !multiplier_context || !work ||
        byte_count > static_cast<uint64_t>(SIZE_MAX) ||
        (gf16 && (byte_count & 1U) != 0))
        return false;
    if (byte_count == 0)
        return true;
    for (uint32_t i = 0; i < transform_size; ++i)
        if (!work[i])
            return false;
    return true;
}

static SparseExecuteContext MakeSparseExecuteContext(
    const leopard::backend::Ops& ops,
    uint64_t byte_count,
    uint16_t zero_multiplier_log,
    const uint8_t* operation_masks,
    size_t operation_count,
    uint32_t skipped_distance,
    PrunedMultiplierLogProvider multiplier_log,
    const void* multiplier_context,
    void** work)
{
    const bool gf16 = zero_multiplier_log == 65535U;
    const PrunedExecutionOps selected = {
        gf16 ? ops.ff16_multiply : ops.ff8_multiply,
        gf16 ? ops.ff16_multiply_add : ops.ff8_multiply_add,
        gf16 ? ops.ff16_fft_butterfly2 : ops.ff8_fft_butterfly2,
        gf16 ? ops.ff16_fft_butterfly4 : ops.ff8_fft_butterfly4,
        gf16 ? ops.ff16_fft_butterfly2_out : ops.ff8_fft_butterfly2_out
    };
    SparseExecuteContext context = {
        &ops,
        selected,
        zero_multiplier_log,
        byte_count,
        operation_masks,
        operation_count,
        0,
        skipped_distance,
        multiplier_log,
        multiplier_context,
        work,
        selected.multiply != NULL && selected.multiply_add != NULL &&
            selected.butterfly != NULL &&
            selected.butterfly_four != NULL &&
            selected.butterfly_out != NULL
    };
    return context;
}

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
    void** work)
{
    if (!ValidateSparseExecution(byte_count, transform_size, shift,
            zero_multiplier_log, operation_masks, retained_bytes,
            multiplier_log, multiplier_context, work))
        return false;
    if (byte_count == 0)
        return true;
    const size_t full_count = SparseForwardButterflyCount(transform_size);
    SparseExecuteContext context = MakeSparseExecuteContext(
        ops, byte_count, zero_multiplier_log, operation_masks, full_count, 0,
        multiplier_log, multiplier_context, work);
    ExecuteSparseNode(context, 0, transform_size, shift);
    return context.valid && context.operation_index == full_count;
}

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
    void** work)
{
    if (!ValidateSparseExecution(byte_count, transform_size, shift,
            zero_multiplier_log, operation_masks, retained_bytes,
            multiplier_log, multiplier_context, work) || !source)
        return false;
    if (byte_count == 0)
        return true;
    for (uint32_t i = 0; i < transform_size; ++i)
        if (!source[i])
            return false;
    const size_t full_count = SparseForwardButterflyCount(transform_size);
    SparseExecuteContext context = MakeSparseExecuteContext(
        ops, byte_count, zero_multiplier_log, operation_masks, full_count,
        transform_size >> 1, multiplier_log, multiplier_context, work);
    if (!ExecuteSparseRootFromSources(
            context, transform_size, shift, source))
        return false;
    ExecuteSparseNode(context, 0, transform_size, shift);
    return context.valid && context.operation_index == full_count;
}

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
    PrunedTransformPlan& plan)
{
    if (!IsPowerOfTwo(field_order) || field_order > 65536U ||
        static_cast<uint32_t>(zero_multiplier_log) + 1U != field_order ||
        !IsPowerOfTwo(size) || size < 2 || size > field_order ||
        shift > field_order - size || (shift & (size - 1)) != 0 ||
        !input_mask || !output_mask || !multiplier_log)
        return false;

    for (uint32_t i = 0; i < size; ++i)
        if (input_mask[i] > 1 || output_mask[i] > 1)
            return false;

    try
    {
        PrunedTransformPlan candidate;
        candidate.size = size;
        candidate.shift = shift;
        candidate.zero_multiplier_log = zero_multiplier_log;
        candidate.inverse = inverse;
        candidate.input_mask.assign(input_mask, input_mask + size);
        candidate.output_mask.assign(output_mask, output_mask + size);

        std::vector<RawPrunedOperation> raw;
        const uint8_t log2_size = Log2PowerOfTwo(size);
        candidate.full_butterfly_count =
            static_cast<size_t>(size >> 1) * log2_size;
        raw.reserve(candidate.full_butterfly_count);
        BuildRawPrunedOperations(
            0, size, shift, inverse,
            multiplier_log, multiplier_context, raw);
        if (raw.size() != candidate.full_butterfly_count)
            return false;
        if (raw.empty())
            return false;
        candidate.first_layer_multiplier_log = multiplier_log(
            multiplier_context, shift + (size >> 1));
        if (candidate.first_layer_multiplier_log > zero_multiplier_log)
            return false;
        for (size_t i = 0; i < raw.size(); ++i)
            if (raw[i].multiplier_log > zero_multiplier_log)
                return false;

        std::vector<uint8_t> live(candidate.input_mask);
        std::vector<uint8_t> live_before(raw.size() * 2U, 0);
        for (size_t i = 0; i < raw.size(); ++i)
        {
            const RawPrunedOperation& operation = raw[i];
            const bool live_x = live[operation.x] != 0;
            const bool live_y = live[operation.y] != 0;
            live_before[i * 2] = static_cast<uint8_t>(live_x);
            live_before[i * 2 + 1] = static_cast<uint8_t>(live_y);
            const bool multiplier_zero =
                operation.multiplier_log == zero_multiplier_log;
            const bool multiplier_one = operation.multiplier_log == 0;
            if (inverse)
            {
                live[operation.x] = static_cast<uint8_t>(
                    (live_x && !multiplier_one) ||
                    (live_y && !multiplier_zero));
                live[operation.y] = static_cast<uint8_t>(live_x || live_y);
            }
            else
            {
                live[operation.x] = static_cast<uint8_t>(
                    live_x || (live_y && !multiplier_zero));
                live[operation.y] = static_cast<uint8_t>(
                    live_x || (live_y && !multiplier_one));
            }
        }

        std::vector<uint8_t> needed(size, 0);
        for (uint32_t i = 0; i < size; ++i)
        {
            needed[i] = static_cast<uint8_t>(
                candidate.output_mask[i] != 0 && live[i] != 0);
            if (candidate.output_mask[i] != 0 && live[i] == 0)
                candidate.zero_outputs.push_back(i);
        }

        std::vector<PrunedTransformOperation> planned(raw.size());
        for (size_t reverse = raw.size(); reverse-- > 0;)
        {
            const RawPrunedOperation& operation = raw[reverse];
            const bool live_x = live_before[reverse * 2] != 0;
            const bool live_y = live_before[reverse * 2 + 1] != 0;
            const bool need_x = needed[operation.x] != 0;
            const bool need_y = needed[operation.y] != 0;
            const bool multiplier_zero =
                operation.multiplier_log == zero_multiplier_log;
            const bool multiplier_one = operation.multiplier_log == 0;
            const bool a_nonzero =
                CoefficientAIsNonzero(inverse, multiplier_one);
            const bool d_nonzero =
                CoefficientDIsNonzero(inverse, multiplier_one);
            const bool a_one =
                CoefficientAIsOne(inverse, multiplier_zero);
            const bool d_one =
                CoefficientDIsOne(inverse, multiplier_zero);

            const bool other_x_term = live_y && !multiplier_zero;
            const bool other_y_term = live_x;
            const bool write_x = need_x && !(
                (live_x && a_one && !other_x_term) ||
                (!live_x && !other_x_term));
            const bool write_y = need_y && !(
                (live_y && d_one && !other_y_term) ||
                (!live_y && !other_y_term));

            uint8_t flags = 0;
            if (live_x) flags |= PrunedLiveX;
            if (live_y) flags |= PrunedLiveY;
            if (need_x) flags |= PrunedNeedX;
            if (need_y) flags |= PrunedNeedY;
            if (write_x) flags |= PrunedWriteX;
            if (write_y) flags |= PrunedWriteY;
            const PrunedTransformOperation result = {
                operation.x,
                operation.y,
                operation.multiplier_log,
                flags
            };
            planned[reverse] = result;

            needed[operation.x] = static_cast<uint8_t>(live_x && (
                (need_x && a_nonzero) || need_y));
            needed[operation.y] = static_cast<uint8_t>(live_y && (
                (need_x && !multiplier_zero) ||
                (need_y && d_nonzero)));
        }

        for (uint32_t i = 0; i < size; ++i)
            if (needed[i] != 0 && candidate.input_mask[i] == 0)
                return false;

        size_t retained_operation_count = 0;
        for (size_t i = 0; i < planned.size(); ++i)
        {
            const uint8_t flags = planned[i].flags;
            if ((flags & (PrunedWriteX | PrunedWriteY)) != 0)
                ++retained_operation_count;
        }
        candidate.operations.reserve(retained_operation_count);
        for (size_t i = 0; i < planned.size(); ++i)
        {
            const PrunedTransformOperation& operation = planned[i];
            const bool write_x = 0 != (operation.flags & PrunedWriteX);
            const bool write_y = 0 != (operation.flags & PrunedWriteY);
            if (!write_x && !write_y)
                continue;
            candidate.operations.push_back(operation);
            const bool need_x = 0 != (operation.flags & PrunedNeedX);
            const bool need_y = 0 != (operation.flags & PrunedNeedY);
            if (need_x != need_y)
                ++candidate.one_output_butterflies;
            const bool live_x = 0 != (operation.flags & PrunedLiveX);
            const bool live_y = 0 != (operation.flags & PrunedLiveY);
            if (live_x != live_y)
                ++candidate.input_zero_specializations;
            if (operation.multiplier_log == zero_multiplier_log)
                ++candidate.zero_multiplier_butterflies;
            if (operation.multiplier_log == 0)
                ++candidate.one_multiplier_butterflies;
        }

        for (size_t i = 0; i < candidate.operations.size();)
        {
            FusedFourMatch match;
            if (MatchFusedFour(candidate.operations, i, inverse, match))
            {
                candidate.fused_four_starts.push_back(
                    static_cast<uint32_t>(i));
                i += 4;
            }
            else
                ++i;
        }

        if (inverse)
        {
            uint32_t prefix = 0;
            while (prefix < size && candidate.input_mask[prefix] != 0)
                ++prefix;
            bool strict_prefix = prefix != 0 && prefix < size;
            for (uint32_t i = prefix; i < size; ++i)
                strict_prefix = strict_prefix && candidate.input_mask[i] == 0;
            for (uint32_t i = 0; i < size; ++i)
                strict_prefix = strict_prefix && candidate.output_mask[i] == 1;
            if (strict_prefix)
            {
                candidate.inverse_source_prefix = prefix;
                candidate.inverse_source_groups.reserve(prefix / 4U);
                for (uint32_t base = 0; base + 4U <= prefix; base += 4U)
                {
                    const PrunedInverseSourceGroup group = {
                        multiplier_log(multiplier_context,
                            shift + base + 1U),
                        multiplier_log(multiplier_context,
                            shift + base + 3U),
                        multiplier_log(multiplier_context,
                            shift + base + 2U)
                    };
                    if (group.multiplier_log01 > zero_multiplier_log ||
                        group.multiplier_log23 > zero_multiplier_log ||
                        group.multiplier_log02 > zero_multiplier_log)
                        return false;
                    candidate.inverse_source_groups.push_back(group);
                }
            }
        }

        plan = std::move(candidate);
        return true;
    }
    catch (...)
    {
        return false;
    }
}

bool ExecutePrunedTransformPlan(
    const leopard::backend::Ops& ops,
    uint64_t byte_count,
    const PrunedTransformPlan& plan,
    void** work)
{
    const bool gf16 = plan.zero_multiplier_log == 65535U;
    if ((plan.zero_multiplier_log != 255U && !gf16) ||
        plan.size < 2 || !IsPowerOfTwo(plan.size) || !work ||
        byte_count > static_cast<uint64_t>(SIZE_MAX) ||
        (gf16 && (byte_count & 1U) != 0))
        return false;

    if (byte_count == 0)
        return true;
    for (uint32_t i = 0; i < plan.size; ++i)
        if (!work[i])
            return false;
    for (size_t i = 0; i < plan.zero_outputs.size(); ++i)
        if (plan.zero_outputs[i] >= plan.size)
            return false;

    const PrunedExecutionOps selected = {
        gf16 ? ops.ff16_multiply : ops.ff8_multiply,
        gf16 ? ops.ff16_multiply_add : ops.ff8_multiply_add,
        plan.inverse
            ? (gf16 ? ops.ff16_ifft_butterfly2 : ops.ff8_ifft_butterfly2)
            : (gf16 ? ops.ff16_fft_butterfly2 : ops.ff8_fft_butterfly2),
        plan.inverse
            ? (gf16 ? ops.ff16_ifft_butterfly4 : ops.ff8_ifft_butterfly4)
            : (gf16 ? ops.ff16_fft_butterfly4 : ops.ff8_fft_butterfly4),
        gf16 ? ops.ff16_fft_butterfly2_out : ops.ff8_fft_butterfly2_out
    };

    size_t fused_index = 0;
    for (size_t i = 0; i < plan.operations.size();)
    {
        if (fused_index < plan.fused_four_starts.size() &&
            plan.fused_four_starts[fused_index] == i)
        {
            FusedFourMatch match;
            if (!MatchFusedFour(plan.operations, i, plan.inverse, match) ||
                match.value0 >= plan.size || match.value1 >= plan.size ||
                match.value2 >= plan.size || match.value3 >= plan.size ||
                match.multiplier_log01 > plan.zero_multiplier_log ||
                match.multiplier_log23 > plan.zero_multiplier_log ||
                match.multiplier_log02 > plan.zero_multiplier_log)
                return false;
            selected.butterfly_four(
                work[match.value0], work[match.value1],
                work[match.value2], work[match.value3],
                match.multiplier_log01, match.multiplier_log23,
                match.multiplier_log02, byte_count);
            ++fused_index;
            i += 4;
            continue;
        }
        const PrunedTransformOperation& operation = plan.operations[i];
        if (operation.x >= plan.size || operation.y >= plan.size ||
            operation.x == operation.y ||
            operation.multiplier_log > plan.zero_multiplier_log ||
            !ExecutePrunedOperation(
                ops, selected, plan.inverse, plan.zero_multiplier_log,
                byte_count, operation, work[operation.x], work[operation.y]))
            return false;
        ++i;
    }
    if (fused_index != plan.fused_four_starts.size())
        return false;
    for (size_t i = 0; i < plan.zero_outputs.size(); ++i)
        memset(work[plan.zero_outputs[i]], 0,
            static_cast<size_t>(byte_count));
    return true;
}

bool ExecutePrunedForwardTransformPlanFromSources(
    const leopard::backend::Ops& ops,
    uint64_t byte_count,
    const PrunedTransformPlan& plan,
    void* const* source,
    void** work)
{
    const bool gf16 = plan.zero_multiplier_log == 65535U;
    if ((plan.zero_multiplier_log != 255U && !gf16) ||
        plan.size < 2 || !IsPowerOfTwo(plan.size) || plan.inverse ||
        !source || !work ||
        plan.input_mask.size() != plan.size ||
        plan.first_layer_multiplier_log > plan.zero_multiplier_log ||
        byte_count > static_cast<uint64_t>(SIZE_MAX) ||
        (gf16 && (byte_count & 1U) != 0))
        return false;

    for (uint32_t i = 0; i < plan.size; ++i)
        if (plan.input_mask[i] != 1)
            return false;
    if (byte_count == 0)
        return true;
    for (uint32_t i = 0; i < plan.size; ++i)
        if (!source[i] || !work[i])
            return false;
    for (size_t i = 0; i < plan.zero_outputs.size(); ++i)
        if (plan.zero_outputs[i] >= plan.size)
            return false;

    const PrunedExecutionOps selected = {
        gf16 ? ops.ff16_multiply : ops.ff8_multiply,
        gf16 ? ops.ff16_multiply_add : ops.ff8_multiply_add,
        gf16 ? ops.ff16_fft_butterfly2 : ops.ff8_fft_butterfly2,
        gf16 ? ops.ff16_fft_butterfly4 : ops.ff8_fft_butterfly4,
        gf16 ? ops.ff16_fft_butterfly2_out : ops.ff8_fft_butterfly2_out
    };
    if (!selected.butterfly_out)
        return false;

    // The forward LCH schedule begins with the root radix-2 layer.  Execute
    // it directly from the immutable coefficient block.  All remaining
    // operations have a smaller coordinate distance and can resume in place.
    const uint32_t half = plan.size >> 1;
    for (uint32_t i = 0; i < half; ++i)
        selected.butterfly_out(
            source[i], source[i + half], work[i], work[i + half],
            plan.first_layer_multiplier_log, byte_count);

    size_t fused_index = 0;
    for (size_t i = 0; i < plan.operations.size();)
    {
        // A complete fused descriptor at the root contains the two operations
        // already executed above plus the following radix-2 layer.  Consume
        // the descriptor metadata, skip only the root operations, and execute
        // its second layer through the ordinary audited operation path.
        if (fused_index < plan.fused_four_starts.size() &&
            plan.fused_four_starts[fused_index] == i &&
            plan.operations[i].x < plan.operations[i].y &&
            plan.operations[i].y - plan.operations[i].x == half)
            ++fused_index;

        const PrunedTransformOperation& operation = plan.operations[i];
        if (operation.x >= plan.size || operation.y >= plan.size ||
            operation.x >= operation.y ||
            operation.multiplier_log > plan.zero_multiplier_log)
            return false;
        if (operation.y - operation.x == half)
        {
            ++i;
            continue;
        }

        if (fused_index < plan.fused_four_starts.size() &&
            plan.fused_four_starts[fused_index] == i)
        {
            FusedFourMatch match;
            if (!MatchFusedFour(plan.operations, i, false, match) ||
                match.value0 >= plan.size || match.value1 >= plan.size ||
                match.value2 >= plan.size || match.value3 >= plan.size ||
                match.multiplier_log01 > plan.zero_multiplier_log ||
                match.multiplier_log23 > plan.zero_multiplier_log ||
                match.multiplier_log02 > plan.zero_multiplier_log)
                return false;
            selected.butterfly_four(
                work[match.value0], work[match.value1],
                work[match.value2], work[match.value3],
                match.multiplier_log01, match.multiplier_log23,
                match.multiplier_log02, byte_count);
            ++fused_index;
            i += 4;
            continue;
        }
        if (!ExecutePrunedOperation(
                ops, selected, false, plan.zero_multiplier_log, byte_count,
                operation, work[operation.x], work[operation.y]))
            return false;
        ++i;
    }
    if (fused_index != plan.fused_four_starts.size())
        return false;
    for (size_t i = 0; i < plan.zero_outputs.size(); ++i)
        memset(work[plan.zero_outputs[i]], 0,
            static_cast<size_t>(byte_count));
    return true;
}

static bool IsConsumedInverseSourceOperation(
    const PrunedTransformOperation& operation,
    uint32_t staged_prefix)
{
    if (operation.x >= staged_prefix || operation.y >= staged_prefix ||
        (operation.x >> 2) != (operation.y >> 2))
        return false;
    const uint32_t x = operation.x & 3U;
    const uint32_t y = operation.y & 3U;
    return (x == 0U && (y == 1U || y == 2U)) ||
        (x == 1U && y == 3U) || (x == 2U && y == 3U);
}

bool ExecutePrunedInverseTransformPlanFromSources(
    const leopard::backend::Ops& ops,
    uint64_t byte_count,
    const PrunedTransformPlan& plan,
    const void* const* source,
    void** work)
{
    const bool gf16 = plan.zero_multiplier_log == 65535U;
    if ((plan.zero_multiplier_log != 255U && !gf16) ||
        plan.size < 2 || !IsPowerOfTwo(plan.size) || !plan.inverse ||
        !source || !work ||
        plan.input_mask.size() != plan.size ||
        plan.output_mask.size() != plan.size ||
        plan.inverse_source_prefix == 0 ||
        plan.inverse_source_prefix >= plan.size ||
        plan.inverse_source_groups.size() !=
            plan.inverse_source_prefix / 4U ||
        byte_count > static_cast<uint64_t>(SIZE_MAX) ||
        (gf16 && (byte_count & 1U) != 0))
        return false;

    for (uint32_t i = 0; i < plan.size; ++i)
    {
        const uint8_t expected =
            static_cast<uint8_t>(i < plan.inverse_source_prefix);
        if (plan.input_mask[i] != expected || plan.output_mask[i] != 1)
            return false;
    }
    for (size_t i = 0; i < plan.inverse_source_groups.size(); ++i)
    {
        const PrunedInverseSourceGroup& group =
            plan.inverse_source_groups[i];
        if (group.multiplier_log01 > plan.zero_multiplier_log ||
            group.multiplier_log23 > plan.zero_multiplier_log ||
            group.multiplier_log02 > plan.zero_multiplier_log)
            return false;
    }
    for (size_t i = 0; i < plan.zero_outputs.size(); ++i)
        if (plan.zero_outputs[i] >= plan.size)
            return false;
    if (byte_count == 0)
        return true;
    for (uint32_t i = 0; i < plan.inverse_source_prefix; ++i)
        if (!source[i])
            return false;
    for (uint32_t i = 0; i < plan.size; ++i)
        if (!work[i])
            return false;

    const PrunedExecutionOps selected = {
        gf16 ? ops.ff16_multiply : ops.ff8_multiply,
        gf16 ? ops.ff16_multiply_add : ops.ff8_multiply_add,
        gf16 ? ops.ff16_ifft_butterfly2 : ops.ff8_ifft_butterfly2,
        gf16 ? ops.ff16_ifft_butterfly4 : ops.ff8_ifft_butterfly4,
        gf16 ? ops.ff16_fft_butterfly2_out : ops.ff8_fft_butterfly2_out
    };
    const leopard::backend::IFFTButterfly4Out butterfly_four_out =
        gf16 ? ops.ff16_ifft_butterfly4_out : ops.ff8_ifft_butterfly4_out;
    if (!selected.multiply || !selected.multiply_add ||
        !selected.butterfly || !selected.butterfly_four ||
        !butterfly_four_out)
        return false;

    const uint32_t staged_prefix = static_cast<uint32_t>(
        plan.inverse_source_groups.size() * 4U);
    for (uint32_t base = 0; base < staged_prefix; base += 4U)
    {
        const PrunedInverseSourceGroup& group =
            plan.inverse_source_groups[base >> 2];
        butterfly_four_out(
            source[base], source[base + 1U],
            source[base + 2U], source[base + 3U],
            work[base], work[base + 1U],
            work[base + 2U], work[base + 3U],
            group.multiplier_log01, group.multiplier_log23,
            group.multiplier_log02, byte_count);
    }
    for (uint32_t i = staged_prefix; i < plan.inverse_source_prefix; ++i)
        memcpy(work[i], source[i], static_cast<size_t>(byte_count));

    size_t fused_index = 0;
    for (size_t i = 0; i < plan.operations.size();)
    {
        while (fused_index < plan.fused_four_starts.size() &&
            plan.fused_four_starts[fused_index] < i)
            ++fused_index;

        const PrunedTransformOperation& operation = plan.operations[i];
        if (operation.x >= plan.size || operation.y >= plan.size ||
            operation.x == operation.y ||
            operation.multiplier_log > plan.zero_multiplier_log)
            return false;
        if (IsConsumedInverseSourceOperation(operation, staged_prefix))
        {
            ++i;
            continue;
        }

        if (fused_index < plan.fused_four_starts.size() &&
            plan.fused_four_starts[fused_index] == i)
        {
            if (plan.operations.size() - i < 4)
                return false;
            bool consumed = false;
            for (size_t j = 0; j < 4; ++j)
                consumed = consumed || IsConsumedInverseSourceOperation(
                    plan.operations[i + j], staged_prefix);
            if (!consumed)
            {
                FusedFourMatch match;
                if (!MatchFusedFour(plan.operations, i, true, match) ||
                    match.value0 >= plan.size || match.value1 >= plan.size ||
                    match.value2 >= plan.size || match.value3 >= plan.size ||
                    match.multiplier_log01 > plan.zero_multiplier_log ||
                    match.multiplier_log23 > plan.zero_multiplier_log ||
                    match.multiplier_log02 > plan.zero_multiplier_log)
                    return false;
                selected.butterfly_four(
                    work[match.value0], work[match.value1],
                    work[match.value2], work[match.value3],
                    match.multiplier_log01, match.multiplier_log23,
                    match.multiplier_log02, byte_count);
                ++fused_index;
                i += 4;
                continue;
            }
        }

        if (!ExecutePrunedOperation(
                ops, selected, true, plan.zero_multiplier_log, byte_count,
                operation, work[operation.x], work[operation.y]))
            return false;
        ++i;
    }
    while (fused_index < plan.fused_four_starts.size() &&
        plan.fused_four_starts[fused_index] < plan.operations.size())
        ++fused_index;
    if (fused_index != plan.fused_four_starts.size())
        return false;
    for (size_t i = 0; i < plan.zero_outputs.size(); ++i)
        memset(work[plan.zero_outputs[i]], 0,
            static_cast<size_t>(byte_count));
    return true;
}

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
