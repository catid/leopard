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
    if (!inverse)
    {
        const uint16_t log_m = multiplier_log(
            multiplier_context, coset + half);
        for (uint32_t offset = 0; offset < half; ++offset)
        {
            const RawPrunedOperation operation = {
                start + offset,
                start + half + offset,
                log_m
            };
            operations.push_back(operation);
        }
    }

    BuildRawPrunedOperations(
        start, half, coset, inverse,
        multiplier_log, multiplier_context, operations);
    BuildRawPrunedOperations(
        start + half, half, coset + half, inverse,
        multiplier_log, multiplier_context, operations);

    if (inverse)
    {
        const uint16_t log_m = multiplier_log(
            multiplier_context, coset + half);
        for (uint32_t offset = 0; offset < half; ++offset)
        {
            const RawPrunedOperation operation = {
                start + offset,
                start + half + offset,
                log_m
            };
            operations.push_back(operation);
        }
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
};

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

} // namespace

namespace leopard2_internal {

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

        candidate.operations.reserve(raw.size());
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
            : (gf16 ? ops.ff16_fft_butterfly2 : ops.ff8_fft_butterfly2)
    };

    for (size_t i = 0; i < plan.operations.size(); ++i)
    {
        const PrunedTransformOperation& operation = plan.operations[i];
        if (operation.x >= plan.size || operation.y >= plan.size ||
            operation.x == operation.y ||
            (!gf16 && operation.multiplier_log > 255U) ||
            !ExecutePrunedOperation(
                ops, selected, plan.inverse, plan.zero_multiplier_log,
                byte_count, operation, work[operation.x], work[operation.y]))
            return false;
    }
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
