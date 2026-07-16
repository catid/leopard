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

#include "leopard2.h"

#include "LeopardCommon.h"
#include "LeopardFF8.h"
#include "LeopardFF16.h"
#include "leopard.h"

#include <algorithm>
#include <limits>
#include <mutex>
#include <new>
#include <string.h>
#include <thread>
#include <vector>

struct leo2_context
{
    leo2_backend backend;
    uint32_t thread_count;
};

struct leo2_codec
{
    leo2_context* context;
    uint32_t original_count;
    uint32_t recovery_count;
    uint32_t parent_count;
    uint32_t padded_side;
    uint32_t parent_dimension;
    leo2_profile profile;
    leo2_field field;
};

struct leo2_decode_plan
{
    const leo2_codec* codec;
    std::vector<uint8_t> original_present;
    std::vector<uint8_t> recovery_present;
    std::vector<uint8_t> coordinate_erased;
    std::vector<uint8_t> requested;
    std::vector<uint8_t> locator8;
    std::vector<uint16_t> locator16;
    uint32_t missing_original_count;
    bool no_op;
    bool direct_xor;
    bool direct_copy;
};

namespace {

static const size_t kScratchAlignment = 64;

struct AddressRange
{
    uintptr_t begin;
    uintptr_t end;
};

struct ScratchLayout
{
    size_t range_offset;
    size_t pointer_offset;
    size_t data_offset;
    size_t total_bytes;
};

static bool CheckedAdd(size_t a, size_t b, size_t& result)
{
    if (a > std::numeric_limits<size_t>::max() - b)
        return false;
    result = a + b;
    return true;
}

static bool CheckedMultiply(size_t a, size_t b, size_t& result)
{
    if (a != 0 && b > std::numeric_limits<size_t>::max() / a)
        return false;
    result = a * b;
    return true;
}

static bool AlignUp(size_t value, size_t alignment, size_t& result)
{
    const size_t mask = alignment - 1;
    size_t sum = 0;
    if (!CheckedAdd(value, mask, sum))
        return false;
    result = sum & ~mask;
    return true;
}

static bool RoundShardBytes(uint64_t shard_bytes, size_t& rounded)
{
    if (shard_bytes == 0 || shard_bytes > std::numeric_limits<size_t>::max())
        return false;
    return AlignUp(static_cast<size_t>(shard_bytes), kScratchAlignment, rounded);
}

static bool ComputeScratchLayout(
    size_t range_count,
    size_t pointer_count,
    size_t slot_count,
    size_t rounded_bytes,
    ScratchLayout& layout)
{
    size_t ranges_bytes = 0;
    size_t pointers_bytes = 0;
    size_t slots_bytes = 0;
    if (!CheckedMultiply(range_count, sizeof(AddressRange), ranges_bytes) ||
        !CheckedMultiply(pointer_count, sizeof(void*), pointers_bytes) ||
        !CheckedMultiply(slot_count, rounded_bytes, slots_bytes))
        return false;

    layout.range_offset = 0;
    if (!AlignUp(ranges_bytes, alignof(void*), layout.pointer_offset))
        return false;
    size_t end_pointers = 0;
    if (!CheckedAdd(layout.pointer_offset, pointers_bytes, end_pointers) ||
        !AlignUp(end_pointers, kScratchAlignment, layout.data_offset) ||
        !CheckedAdd(layout.data_offset, slots_bytes, layout.total_bytes))
        return false;
    return true;
}

static uint32_t CeilPow2(uint64_t value)
{
    if (value == 0 || value > 65536)
        return 0;
    uint32_t result = 1;
    while (result < value)
        result <<= 1;
    return result;
}

static uint32_t CoordinateForOriginal(const leo2_codec* codec, uint32_t index)
{
    return codec->profile == LEO2_PROFILE_LEGACY_HIGH_V1
        ? codec->padded_side + index
        : index;
}

static uint32_t CoordinateForRecovery(const leo2_codec* codec, uint32_t index)
{
    return codec->profile == LEO2_PROFILE_LEGACY_HIGH_V1
        ? index
        : codec->padded_side + index;
}

static leo2_backend RuntimeBackend()
{
#if defined(LEO_TRY_NEON)
    if (leopard::CpuHasNeon)
        return LEO2_BACKEND_NEON;
#endif
#if defined(LEO_TRY_AVX2)
    if (leopard::CpuHasAVX2)
        return LEO2_BACKEND_AVX2;
#endif
#if !defined(LEO_TARGET_MOBILE) || defined(LEO_USE_SSE2NEON)
    if (leopard::CpuHasSSSE3)
        return LEO2_BACKEND_SSSE3;
#endif
    return LEO2_BACKEND_SCALAR;
}

static leo2_result EnsureInitialized()
{
    static std::mutex mutex;
    static bool attempted = false;
    static int result = Leopard_CallInitialize;
    std::lock_guard<std::mutex> lock(mutex);
    if (!attempted)
    {
        result = leo_init();
        attempted = true;
    }
    return result == Leopard_Success ? LEO2_SUCCESS : LEO2_INTERNAL_ERROR;
}

static bool MakeRange(const void* pointer, uint64_t bytes, AddressRange& range)
{
    if (!pointer || bytes == 0 || bytes > std::numeric_limits<uintptr_t>::max())
        return false;
    const uintptr_t begin = reinterpret_cast<uintptr_t>(pointer);
    const uintptr_t length = static_cast<uintptr_t>(bytes);
    if (begin > std::numeric_limits<uintptr_t>::max() - length)
        return false;
    range.begin = begin;
    range.end = begin + length;
    return true;
}

static bool RangesOverlap(const AddressRange& a, const AddressRange& b)
{
    return a.begin < b.end && b.begin < a.end;
}

static bool RangeLess(const AddressRange& a, const AddressRange& b)
{
    if (a.begin != b.begin)
        return a.begin < b.begin;
    return a.end < b.end;
}

static size_t MergeRanges(AddressRange* ranges, size_t count)
{
    if (count == 0)
        return 0;
    std::sort(ranges, ranges + count, RangeLess);
    size_t merged = 1;
    for (size_t i = 1; i < count; ++i)
    {
        AddressRange& back = ranges[merged - 1];
        if (ranges[i].begin <= back.end)
        {
            if (ranges[i].end > back.end)
                back.end = ranges[i].end;
        }
        else
            ranges[merged++] = ranges[i];
    }
    return merged;
}

static leo2_result ValidateDisjointRanges(
    AddressRange* input_ranges,
    size_t input_count,
    AddressRange* output_ranges,
    size_t output_count)
{
    input_count = MergeRanges(input_ranges, input_count);
    std::sort(output_ranges, output_ranges + output_count, RangeLess);
    for (size_t i = 1; i < output_count; ++i)
        if (RangesOverlap(output_ranges[i - 1], output_ranges[i]))
            return LEO2_OVERLAP;

    size_t input_i = 0;
    for (size_t output_i = 0; output_i < output_count; ++output_i)
    {
        while (input_i < input_count && input_ranges[input_i].end <= output_ranges[output_i].begin)
            ++input_i;
        if (input_i < input_count && RangesOverlap(input_ranges[input_i], output_ranges[output_i]))
            return LEO2_OVERLAP;
    }
    return LEO2_SUCCESS;
}

static leo2_result CheckScratch(
    void* scratch,
    size_t scratch_bytes,
    const ScratchLayout& layout,
    AddressRange& scratch_range)
{
    if (scratch_bytes < layout.total_bytes)
        return LEO2_SCRATCH_TOO_SMALL;
    if (layout.total_bytes == 0)
        return LEO2_SUCCESS;
    if (!scratch)
        return LEO2_INVALID_ARGUMENT;
    if ((reinterpret_cast<uintptr_t>(scratch) & (kScratchAlignment - 1)) != 0)
        return LEO2_BAD_ALIGNMENT;
    if (!MakeRange(scratch, layout.total_bytes, scratch_range))
        return LEO2_INVALID_ARGUMENT;
    return LEO2_SUCCESS;
}

static leo2_result ValidateEncodeBuffers(
    const leo2_codec* codec,
    uint64_t shard_bytes,
    const void* const* original,
    void* const* recovery,
    void* scratch,
    const ScratchLayout& layout)
{
    if (!original || !recovery)
        return LEO2_INVALID_ARGUMENT;
    AddressRange scratch_range;
    if (!MakeRange(scratch, layout.total_bytes, scratch_range))
        return LEO2_INVALID_ARGUMENT;

    for (uint32_t i = 0; i < codec->original_count; ++i)
    {
        AddressRange range;
        if (!MakeRange(original[i], shard_bytes, range))
            return LEO2_INVALID_ARGUMENT;
        if (RangesOverlap(range, scratch_range))
            return LEO2_OVERLAP;
    }
    for (uint32_t i = 0; i < codec->recovery_count; ++i)
    {
        if (!recovery[i])
            continue;
        AddressRange range;
        if (!MakeRange(recovery[i], shard_bytes, range))
            return LEO2_INVALID_ARGUMENT;
        if (RangesOverlap(range, scratch_range))
            return LEO2_OVERLAP;
    }

    AddressRange* ranges = reinterpret_cast<AddressRange*>(scratch);
    size_t input_count = 0;
    size_t output_count = 0;
    for (uint32_t i = 0; i < codec->original_count; ++i)
        MakeRange(original[i], shard_bytes, ranges[input_count++]);
    AddressRange* outputs = ranges + codec->original_count;
    for (uint32_t i = 0; i < codec->recovery_count; ++i)
        if (recovery[i])
            MakeRange(recovery[i], shard_bytes, outputs[output_count++]);
    return ValidateDisjointRanges(ranges, input_count, outputs, output_count);
}

static leo2_result EncodeLayout(
    const leo2_codec* codec,
    uint64_t shard_bytes,
    ScratchLayout& layout,
    size_t& rounded_bytes,
    size_t& work_count)
{
    if (!codec || !RoundShardBytes(shard_bytes, rounded_bytes))
        return LEO2_INVALID_ARGUMENT;
    /*
        Legacy GF16 stores 32 symbols per 64-byte ALTMAP tile.  Truncating a
        tile can discard a nonzero high half that is required for decoding, so
        partial tiles need a separately versioned compact-tail construction.
    */
    if (codec->field == LEO2_FIELD_GF16 && (shard_bytes & 63u) != 0)
        return LEO2_UNSUPPORTED;
    work_count = static_cast<size_t>(codec->padded_side) * 2;
    const size_t parity_slots = codec->profile == LEO2_PROFILE_LOW_V1
        ? codec->recovery_count : 0;
    const size_t range_count = static_cast<size_t>(codec->original_count) + codec->recovery_count;
    const size_t pointer_count = static_cast<size_t>(codec->original_count) + work_count + parity_slots;
    const size_t slot_count = pointer_count;
    if (!ComputeScratchLayout(range_count, pointer_count, slot_count, rounded_bytes, layout))
        return LEO2_INVALID_COUNTS;
    return LEO2_SUCCESS;
}

static leo2_result DecodeLayout(
    const leo2_codec* codec,
    uint64_t shard_bytes,
    ScratchLayout& layout,
    size_t& rounded_bytes)
{
    if (!codec || !RoundShardBytes(shard_bytes, rounded_bytes))
        return LEO2_INVALID_ARGUMENT;
    if (codec->field == LEO2_FIELD_GF16 && (shard_bytes & 63u) != 0)
        return LEO2_UNSUPPORTED;
    const size_t range_count = static_cast<size_t>(codec->original_count) * 2 + codec->recovery_count;
    const size_t pointer_count = static_cast<size_t>(codec->parent_count) * 2;
    const size_t slot_count = static_cast<size_t>(codec->original_count) +
        codec->recovery_count + codec->parent_count;
    if (!ComputeScratchLayout(range_count, pointer_count, slot_count, rounded_bytes, layout))
        return LEO2_INVALID_COUNTS;
    return LEO2_SUCCESS;
}

static void CopyAndPad(void* destination, const void* source, size_t bytes, size_t rounded)
{
    memcpy(destination, source, bytes);
    if (rounded > bytes)
        memset(static_cast<uint8_t*>(destination) + bytes, 0, rounded - bytes);
}

static leo2_result ValidateDecodeBuffers(
    const leo2_decode_plan* plan,
    uint64_t shard_bytes,
    const void* const* original,
    const void* const* recovery,
    void* const* restored,
    void* scratch,
    const ScratchLayout& layout)
{
    const leo2_codec* codec = plan->codec;
    if (!original || !recovery || !restored)
        return LEO2_INVALID_ARGUMENT;
    AddressRange scratch_range;
    if (!MakeRange(scratch, layout.total_bytes, scratch_range))
        return LEO2_INVALID_ARGUMENT;

    size_t input_count = 0;
    size_t output_count = 0;
    for (uint32_t i = 0; i < codec->original_count; ++i)
    {
        if ((original[i] != NULL) != (plan->original_present[i] != 0))
            return LEO2_INVALID_ARGUMENT;
        if (original[i])
        {
            AddressRange range;
            if (!MakeRange(original[i], shard_bytes, range))
                return LEO2_INVALID_ARGUMENT;
            if (RangesOverlap(range, scratch_range))
                return LEO2_OVERLAP;
            ++input_count;
        }
        else
        {
            if (!restored[i])
                return LEO2_INVALID_ARGUMENT;
            AddressRange range;
            if (!MakeRange(restored[i], shard_bytes, range))
                return LEO2_INVALID_ARGUMENT;
            if (RangesOverlap(range, scratch_range))
                return LEO2_OVERLAP;
            ++output_count;
        }
    }
    for (uint32_t i = 0; i < codec->recovery_count; ++i)
    {
        if ((recovery[i] != NULL) != (plan->recovery_present[i] != 0))
            return LEO2_INVALID_ARGUMENT;
        if (recovery[i])
        {
            AddressRange range;
            if (!MakeRange(recovery[i], shard_bytes, range))
                return LEO2_INVALID_ARGUMENT;
            if (RangesOverlap(range, scratch_range))
                return LEO2_OVERLAP;
            ++input_count;
        }
    }

    AddressRange* ranges = reinterpret_cast<AddressRange*>(scratch);
    input_count = 0;
    output_count = 0;
    for (uint32_t i = 0; i < codec->original_count; ++i)
        if (original[i])
            MakeRange(original[i], shard_bytes, ranges[input_count++]);
    for (uint32_t i = 0; i < codec->recovery_count; ++i)
        if (recovery[i])
            MakeRange(recovery[i], shard_bytes, ranges[input_count++]);
    AddressRange* outputs = ranges + input_count;
    for (uint32_t i = 0; i < codec->original_count; ++i)
        if (!original[i])
            MakeRange(restored[i], shard_bytes, outputs[output_count++]);
    return ValidateDisjointRanges(ranges, input_count, outputs, output_count);
}

} // namespace

extern "C" {

LEO2_EXPORT const char* leo2_result_string(leo2_result result)
{
    switch (result)
    {
    case LEO2_SUCCESS: return "Operation succeeded";
    case LEO2_NEED_MORE_DATA: return "Not enough received shards";
    case LEO2_INVALID_ARGUMENT: return "Invalid argument";
    case LEO2_INVALID_COUNTS: return "Invalid or unsupported shard counts";
    case LEO2_UNSUPPORTED: return "Requested profile, field, or backend is unsupported";
    case LEO2_SCRATCH_TOO_SMALL: return "Scratch buffer is too small";
    case LEO2_BAD_ALIGNMENT: return "Scratch buffer has insufficient alignment";
    case LEO2_OVERLAP: return "Unsupported input, output, or scratch overlap";
    case LEO2_OUT_OF_MEMORY: return "Allocation failed during setup";
    case LEO2_INTERNAL_ERROR: return "Internal initialization or execution error";
    }
    return "Unknown Leopard2 result";
}

LEO2_EXPORT leo2_result leo2_context_create(
    const leo2_context_options* options,
    leo2_context** context_out)
{
    if (!context_out)
        return LEO2_INVALID_ARGUMENT;
    *context_out = NULL;
    if (options && options->struct_size < sizeof(leo2_context_options))
        return LEO2_INVALID_ARGUMENT;
    const leo2_result initialized = EnsureInitialized();
    if (initialized != LEO2_SUCCESS)
        return initialized;

    const leo2_backend actual = RuntimeBackend();
    const leo2_backend requested = options ? options->backend : LEO2_BACKEND_AUTO;
    if (requested < LEO2_BACKEND_AUTO || requested > LEO2_BACKEND_NEON)
        return LEO2_INVALID_ARGUMENT;
    if (requested != LEO2_BACKEND_AUTO && requested != actual)
        return LEO2_UNSUPPORTED;

    leo2_context* context = new (std::nothrow) leo2_context;
    if (!context)
        return LEO2_OUT_OF_MEMORY;
    context->backend = actual;
    uint32_t threads = options ? options->thread_count : 0;
    if (threads == 0)
    {
        threads = static_cast<uint32_t>(std::thread::hardware_concurrency());
        if (threads == 0)
            threads = 1;
    }
    context->thread_count = threads;
    *context_out = context;
    return LEO2_SUCCESS;
}

LEO2_EXPORT void leo2_context_destroy(leo2_context* context)
{
    delete context;
}

LEO2_EXPORT leo2_backend leo2_context_backend(const leo2_context* context)
{
    return context ? context->backend : LEO2_BACKEND_AUTO;
}

LEO2_EXPORT uint32_t leo2_context_thread_count(const leo2_context* context)
{
    return context ? context->thread_count : 0;
}

LEO2_EXPORT leo2_result leo2_codec_create(
    leo2_context* context,
    uint32_t original_count,
    uint32_t recovery_count,
    leo2_profile profile,
    leo2_field field,
    const leo2_codec_options* options,
    leo2_codec** codec_out)
{
    if (!context || !codec_out || original_count == 0 || recovery_count == 0)
        return LEO2_INVALID_ARGUMENT;
    *codec_out = NULL;
    if (options)
    {
        if (options->struct_size < sizeof(leo2_codec_options) || options->flags != 0)
            return LEO2_INVALID_ARGUMENT;
    }
    if (profile == LEO2_PROFILE_AUTO)
        profile = recovery_count <= original_count
            ? LEO2_PROFILE_LEGACY_HIGH_V1 : LEO2_PROFILE_LOW_V1;
    if (profile != LEO2_PROFILE_LEGACY_HIGH_V1 && profile != LEO2_PROFILE_LOW_V1)
        return profile == LEO2_PROFILE_EXACT_EXPERIMENTAL_V1
            ? LEO2_UNSUPPORTED : LEO2_INVALID_ARGUMENT;

    const uint32_t padded = CeilPow2(
        profile == LEO2_PROFILE_LEGACY_HIGH_V1 ? recovery_count : original_count);
    if (padded == 0)
        return LEO2_INVALID_COUNTS;
    const uint32_t parent = CeilPow2(static_cast<uint64_t>(padded) +
        (profile == LEO2_PROFILE_LEGACY_HIGH_V1 ? original_count : recovery_count));
    if (parent == 0)
        return LEO2_INVALID_COUNTS;

    if (field == LEO2_FIELD_AUTO)
        field = parent <= leopard::ff8::kOrder ? LEO2_FIELD_GF8 : LEO2_FIELD_GF16;
    if (field == LEO2_FIELD_GF8 && parent > leopard::ff8::kOrder)
        return LEO2_INVALID_COUNTS;
    if (field == LEO2_FIELD_GF16 && parent > leopard::ff16::kOrder)
        return LEO2_INVALID_COUNTS;
    if (field != LEO2_FIELD_GF8 && field != LEO2_FIELD_GF16)
        return LEO2_INVALID_ARGUMENT;

    leo2_codec* codec = new (std::nothrow) leo2_codec;
    if (!codec)
        return LEO2_OUT_OF_MEMORY;
    codec->context = context;
    codec->original_count = original_count;
    codec->recovery_count = recovery_count;
    codec->parent_count = parent;
    codec->padded_side = padded;
    codec->parent_dimension = profile == LEO2_PROFILE_LEGACY_HIGH_V1
        ? parent - padded : padded;
    codec->profile = profile;
    codec->field = field;
    *codec_out = codec;
    return LEO2_SUCCESS;
}

LEO2_EXPORT void leo2_codec_destroy(leo2_codec* codec)
{
    delete codec;
}

LEO2_EXPORT uint32_t leo2_codec_original_count(const leo2_codec* codec)
{
    return codec ? codec->original_count : 0;
}

LEO2_EXPORT uint32_t leo2_codec_recovery_count(const leo2_codec* codec)
{
    return codec ? codec->recovery_count : 0;
}

LEO2_EXPORT uint32_t leo2_codec_parent_count(const leo2_codec* codec)
{
    return codec ? codec->parent_count : 0;
}

LEO2_EXPORT uint32_t leo2_codec_padded_side(const leo2_codec* codec)
{
    return codec ? codec->padded_side : 0;
}

LEO2_EXPORT leo2_profile leo2_codec_profile(const leo2_codec* codec)
{
    return codec ? codec->profile : LEO2_PROFILE_AUTO;
}

LEO2_EXPORT leo2_field leo2_codec_field(const leo2_codec* codec)
{
    return codec ? codec->field : LEO2_FIELD_AUTO;
}

LEO2_EXPORT size_t leo2_scratch_alignment(void)
{
    return kScratchAlignment;
}

LEO2_EXPORT leo2_result leo2_encode_scratch_size(
    const leo2_codec* codec,
    uint64_t shard_bytes,
    size_t* scratch_bytes_out)
{
    if (!scratch_bytes_out)
        return LEO2_INVALID_ARGUMENT;
    ScratchLayout layout;
    size_t rounded = 0;
    size_t work_count = 0;
    const leo2_result result = EncodeLayout(codec, shard_bytes, layout, rounded, work_count);
    if (result != LEO2_SUCCESS)
        return result;
    *scratch_bytes_out = layout.total_bytes;
    return LEO2_SUCCESS;
}

LEO2_EXPORT leo2_result leo2_encode(
    const leo2_codec* codec,
    uint64_t shard_bytes,
    const void* const* original,
    void* const* recovery,
    void* scratch,
    size_t scratch_bytes)
{
    ScratchLayout layout;
    size_t rounded = 0;
    size_t work_count = 0;
    leo2_result result = EncodeLayout(codec, shard_bytes, layout, rounded, work_count);
    if (result != LEO2_SUCCESS)
        return result;
    AddressRange scratch_range;
    result = CheckScratch(scratch, scratch_bytes, layout, scratch_range);
    if (result != LEO2_SUCCESS)
        return result;
    result = ValidateEncodeBuffers(codec, shard_bytes, original, recovery, scratch, layout);
    if (result != LEO2_SUCCESS)
        return result;

    uint8_t* base = static_cast<uint8_t*>(scratch);
    void** pointers = reinterpret_cast<void**>(base + layout.pointer_offset);
    uint8_t* slots = base + layout.data_offset;
    for (uint32_t i = 0; i < codec->original_count; ++i)
    {
        pointers[i] = slots + static_cast<size_t>(i) * rounded;
        CopyAndPad(pointers[i], original[i], static_cast<size_t>(shard_bytes), rounded);
    }
    void** work = pointers + codec->original_count;
    for (size_t i = 0; i < work_count; ++i)
        work[i] = slots + (static_cast<size_t>(codec->original_count) + i) * rounded;

    const void* const* padded_original = const_cast<const void* const*>(pointers);
    if (codec->profile == LEO2_PROFILE_LEGACY_HIGH_V1)
    {
        if (codec->padded_side == 1)
        {
            uint8_t* parity = static_cast<uint8_t*>(work[0]);
            memcpy(parity, pointers[0], rounded);
            for (uint32_t i = 1; i < codec->original_count; ++i)
            {
                const uint8_t* source = static_cast<const uint8_t*>(pointers[i]);
                for (size_t j = 0; j < rounded; ++j)
                    parity[j] ^= source[j];
            }
        }
        else if (codec->field == LEO2_FIELD_GF8)
            leopard::ff8::ReedSolomonEncode(
                rounded, codec->original_count, codec->recovery_count,
                codec->padded_side, padded_original, work);
        else
            leopard::ff16::ReedSolomonEncode(
                rounded, codec->original_count, codec->recovery_count,
                codec->padded_side, padded_original, work);

        for (uint32_t i = 0; i < codec->recovery_count; ++i)
            if (recovery[i])
                memcpy(recovery[i], work[i], static_cast<size_t>(shard_bytes));
    }
    else
    {
        void** parity = work + work_count;
        for (uint32_t i = 0; i < codec->recovery_count; ++i)
        {
            parity[i] = recovery[i]
                ? slots + (static_cast<size_t>(codec->original_count) +
                    work_count + i) * rounded
                : NULL;
        }
        if (codec->padded_side == 1)
        {
            for (uint32_t i = 0; i < codec->recovery_count; ++i)
                if (recovery[i])
                    memcpy(recovery[i], original[0], static_cast<size_t>(shard_bytes));
            return LEO2_SUCCESS;
        }
        if (codec->field == LEO2_FIELD_GF8)
            leopard::ff8::ReedSolomonEncodeLow(
                rounded, codec->original_count, codec->recovery_count,
                codec->padded_side, padded_original, parity, work);
        else
            leopard::ff16::ReedSolomonEncodeLow(
                rounded, codec->original_count, codec->recovery_count,
                codec->padded_side, padded_original, parity, work);
        for (uint32_t i = 0; i < codec->recovery_count; ++i)
            if (recovery[i])
                memcpy(recovery[i], parity[i], static_cast<size_t>(shard_bytes));
    }
    return LEO2_SUCCESS;
}

LEO2_EXPORT leo2_result leo2_encode_batch(
    const leo2_codec* codec,
    const leo2_encode_batch_item* items,
    size_t item_count)
{
    if (item_count != 0 && !items)
        return LEO2_INVALID_ARGUMENT;
    for (size_t i = 0; i < item_count; ++i)
    {
        const leo2_result result = leo2_encode(
            codec, items[i].shard_bytes, items[i].original, items[i].recovery,
            items[i].scratch, items[i].scratch_bytes);
        if (result != LEO2_SUCCESS)
            return result;
    }
    return LEO2_SUCCESS;
}

LEO2_EXPORT leo2_result leo2_decode_plan_create(
    const leo2_codec* codec,
    const uint8_t* original_present,
    const uint8_t* recovery_present,
    leo2_decode_plan** plan_out)
{
    if (!codec || !original_present || !recovery_present || !plan_out)
        return LEO2_INVALID_ARGUMENT;
    *plan_out = NULL;
    uint32_t present_count = 0;
    uint32_t missing_original_count = 0;
    for (uint32_t i = 0; i < codec->original_count; ++i)
    {
        if (original_present[i] > 1)
            return LEO2_INVALID_ARGUMENT;
        if (original_present[i])
            ++present_count;
        else
            ++missing_original_count;
    }
    for (uint32_t i = 0; i < codec->recovery_count; ++i)
    {
        if (recovery_present[i] > 1)
            return LEO2_INVALID_ARGUMENT;
        if (recovery_present[i])
            ++present_count;
    }
    if (present_count < codec->original_count)
        return LEO2_NEED_MORE_DATA;

    leo2_decode_plan* plan = new (std::nothrow) leo2_decode_plan;
    if (!plan)
        return LEO2_OUT_OF_MEMORY;
    try
    {
        plan->codec = codec;
        plan->original_present.assign(original_present, original_present + codec->original_count);
        plan->recovery_present.assign(recovery_present, recovery_present + codec->recovery_count);
        plan->coordinate_erased.assign(codec->parent_count, 0);
        plan->requested.assign(codec->parent_count, 0);
        plan->missing_original_count = missing_original_count;
        plan->no_op = missing_original_count == 0;
        plan->direct_xor = codec->profile == LEO2_PROFILE_LEGACY_HIGH_V1 &&
            codec->padded_side == 1 && missing_original_count == 1;
        plan->direct_copy = codec->profile == LEO2_PROFILE_LOW_V1 &&
            codec->padded_side == 1 && missing_original_count == 1;

        if (codec->profile == LEO2_PROFILE_LEGACY_HIGH_V1)
        {
            for (uint32_t i = codec->recovery_count; i < codec->padded_side; ++i)
                plan->coordinate_erased[i] = 1;
        }
        else
        {
            for (uint32_t i = codec->padded_side + codec->recovery_count;
                 i < codec->parent_count; ++i)
                plan->coordinate_erased[i] = 1;
        }
        for (uint32_t i = 0; i < codec->original_count; ++i)
        {
            if (!original_present[i])
            {
                const uint32_t coordinate = CoordinateForOriginal(codec, i);
                plan->coordinate_erased[coordinate] = 1;
                plan->requested[coordinate] = 1;
            }
        }
        for (uint32_t i = 0; i < codec->recovery_count; ++i)
            if (!recovery_present[i])
                plan->coordinate_erased[CoordinateForRecovery(codec, i)] = 1;

        if (!plan->no_op && !plan->direct_xor && !plan->direct_copy)
        {
            if (codec->field == LEO2_FIELD_GF8)
            {
                plan->locator8.resize(codec->parent_count);
                leopard::ff8::PrepareDecode(codec->parent_count,
                    &plan->coordinate_erased[0], &plan->locator8[0]);
            }
            else
            {
                plan->locator16.resize(codec->parent_count);
                leopard::ff16::PrepareDecode(codec->parent_count,
                    &plan->coordinate_erased[0], &plan->locator16[0]);
            }
        }
    }
    catch (const std::bad_alloc&)
    {
        delete plan;
        return LEO2_OUT_OF_MEMORY;
    }
    *plan_out = plan;
    return LEO2_SUCCESS;
}

LEO2_EXPORT void leo2_decode_plan_destroy(leo2_decode_plan* plan)
{
    delete plan;
}

LEO2_EXPORT uint32_t leo2_decode_plan_missing_original_count(
    const leo2_decode_plan* plan)
{
    return plan ? plan->missing_original_count : 0;
}

LEO2_EXPORT leo2_result leo2_decode_plan_scratch_size(
    const leo2_decode_plan* plan,
    uint64_t shard_bytes,
    size_t* scratch_bytes_out)
{
    if (!plan || !scratch_bytes_out || shard_bytes == 0)
        return LEO2_INVALID_ARGUMENT;
    if (plan->no_op)
    {
        *scratch_bytes_out = 0;
        return LEO2_SUCCESS;
    }
    ScratchLayout layout;
    size_t rounded = 0;
    const leo2_result result = DecodeLayout(plan->codec, shard_bytes, layout, rounded);
    if (result != LEO2_SUCCESS)
        return result;
    *scratch_bytes_out = layout.total_bytes;
    return LEO2_SUCCESS;
}

LEO2_EXPORT leo2_result leo2_decode_plan_execute(
    const leo2_decode_plan* plan,
    uint64_t shard_bytes,
    const void* const* original,
    const void* const* recovery,
    void* const* restored_original,
    void* scratch,
    size_t scratch_bytes)
{
    if (!plan || shard_bytes == 0)
        return LEO2_INVALID_ARGUMENT;
    if (plan->no_op)
        return LEO2_SUCCESS;

    ScratchLayout layout;
    size_t rounded = 0;
    leo2_result result = DecodeLayout(plan->codec, shard_bytes, layout, rounded);
    if (result != LEO2_SUCCESS)
        return result;
    AddressRange scratch_range;
    result = CheckScratch(scratch, scratch_bytes, layout, scratch_range);
    if (result != LEO2_SUCCESS)
        return result;
    result = ValidateDecodeBuffers(
        plan, shard_bytes, original, recovery, restored_original, scratch, layout);
    if (result != LEO2_SUCCESS)
        return result;

    const leo2_codec* codec = plan->codec;
    if (plan->direct_copy)
    {
        for (uint32_t i = 0; i < codec->recovery_count; ++i)
            if (recovery[i])
            {
                memcpy(restored_original[0], recovery[i], static_cast<size_t>(shard_bytes));
                return LEO2_SUCCESS;
            }
        return LEO2_NEED_MORE_DATA;
    }
    if (plan->direct_xor)
    {
        uint32_t missing = 0;
        while (missing < codec->original_count && original[missing])
            ++missing;
        if (missing == codec->original_count || !recovery[0])
            return LEO2_NEED_MORE_DATA;
        uint8_t* output = static_cast<uint8_t*>(restored_original[missing]);
        memcpy(output, recovery[0], static_cast<size_t>(shard_bytes));
        for (uint32_t i = 0; i < codec->original_count; ++i)
        {
            if (!original[i])
                continue;
            const uint8_t* source = static_cast<const uint8_t*>(original[i]);
            for (size_t j = 0; j < static_cast<size_t>(shard_bytes); ++j)
                output[j] ^= source[j];
        }
        return LEO2_SUCCESS;
    }

    uint8_t* base = static_cast<uint8_t*>(scratch);
    void** coordinate_data = reinterpret_cast<void**>(base + layout.pointer_offset);
    void** work = coordinate_data + codec->parent_count;
    uint8_t* slots = base + layout.data_offset;
    std::fill(coordinate_data, coordinate_data + codec->parent_count, static_cast<void*>(NULL));

    for (uint32_t i = 0; i < codec->original_count; ++i)
    {
        uint8_t* slot = slots + static_cast<size_t>(i) * rounded;
        if (original[i])
        {
            CopyAndPad(slot, original[i], static_cast<size_t>(shard_bytes), rounded);
            coordinate_data[CoordinateForOriginal(codec, i)] = slot;
        }
    }
    for (uint32_t i = 0; i < codec->recovery_count; ++i)
    {
        uint8_t* slot = slots + (static_cast<size_t>(codec->original_count) + i) * rounded;
        if (recovery[i])
        {
            CopyAndPad(slot, recovery[i], static_cast<size_t>(shard_bytes), rounded);
            coordinate_data[CoordinateForRecovery(codec, i)] = slot;
        }
    }
    const size_t work_base = static_cast<size_t>(codec->original_count) + codec->recovery_count;
    for (uint32_t i = 0; i < codec->parent_count; ++i)
        work[i] = slots + (work_base + i) * rounded;

    const void* const* coordinate_input = const_cast<const void* const*>(coordinate_data);
    if (codec->field == LEO2_FIELD_GF8)
        leopard::ff8::ReedSolomonDecodePrepared(
            rounded, codec->parent_count, coordinate_input, &plan->requested[0],
            &plan->locator8[0], work);
    else
        leopard::ff16::ReedSolomonDecodePrepared(
            rounded, codec->parent_count, coordinate_input, &plan->requested[0],
            &plan->locator16[0], work);

    for (uint32_t i = 0; i < codec->original_count; ++i)
        if (!original[i])
            memcpy(restored_original[i], work[CoordinateForOriginal(codec, i)],
                static_cast<size_t>(shard_bytes));
    return LEO2_SUCCESS;
}

LEO2_EXPORT leo2_result leo2_decode_plan_execute_batch(
    const leo2_decode_plan* plan,
    const leo2_decode_batch_item* items,
    size_t item_count)
{
    if (item_count != 0 && !items)
        return LEO2_INVALID_ARGUMENT;
    for (size_t i = 0; i < item_count; ++i)
    {
        const leo2_result result = leo2_decode_plan_execute(
            plan, items[i].shard_bytes, items[i].original, items[i].recovery,
            items[i].restored_original, items[i].scratch, items[i].scratch_bytes);
        if (result != LEO2_SUCCESS)
            return result;
    }
    return LEO2_SUCCESS;
}

LEO2_EXPORT leo2_result leo2_decode_scratch_size(
    const leo2_codec* codec,
    uint64_t shard_bytes,
    size_t* scratch_bytes_out)
{
    if (!scratch_bytes_out)
        return LEO2_INVALID_ARGUMENT;
    ScratchLayout layout;
    size_t rounded = 0;
    const leo2_result result = DecodeLayout(codec, shard_bytes, layout, rounded);
    if (result != LEO2_SUCCESS)
        return result;
    *scratch_bytes_out = layout.total_bytes;
    return LEO2_SUCCESS;
}

LEO2_EXPORT leo2_result leo2_decode(
    const leo2_codec* codec,
    uint64_t shard_bytes,
    const uint8_t* original_present,
    const uint8_t* recovery_present,
    const void* const* original,
    const void* const* recovery,
    void* const* restored_original,
    void* scratch,
    size_t scratch_bytes)
{
    leo2_decode_plan* plan = NULL;
    leo2_result result = leo2_decode_plan_create(
        codec, original_present, recovery_present, &plan);
    if (result != LEO2_SUCCESS)
        return result;
    result = leo2_decode_plan_execute(
        plan, shard_bytes, original, recovery, restored_original, scratch, scratch_bytes);
    leo2_decode_plan_destroy(plan);
    return result;
}

} // extern "C"
