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

#include "Leopard2Direct.h"
#include "leopard2.h"

#include <algorithm>
#include <array>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <limits>
#include <new>
#include <stdexcept>
#include <vector>

#if defined(_MSC_VER)
#include <malloc.h>
#endif

namespace {

class AlignedBuffer
{
public:
    explicit AlignedBuffer(size_t bytes)
        : data_(NULL)
        , bytes_(bytes)
    {
        if (bytes_ == 0)
            return;
#if defined(_MSC_VER)
        data_ = _aligned_malloc(bytes_, leo2_scratch_alignment());
#else
        if (posix_memalign(&data_, leo2_scratch_alignment(), bytes_) != 0)
            data_ = NULL;
#endif
        if (!data_)
            throw std::bad_alloc();
        memset(data_, 0, bytes_);
    }

    ~AlignedBuffer()
    {
#if defined(_MSC_VER)
        _aligned_free(data_);
#else
        free(data_);
#endif
    }

    void* data() { return data_; }
    size_t size() const { return bytes_; }

private:
    AlignedBuffer(const AlignedBuffer&);
    AlignedBuffer& operator=(const AlignedBuffer&);

    void* data_;
    size_t bytes_;
};

class CodecOwner
{
public:
    CodecOwner() : codec(NULL) {}
    ~CodecOwner() { leo2_codec_destroy(codec); }
    leo2_codec* codec;

private:
    CodecOwner(const CodecOwner&);
    CodecOwner& operator=(const CodecOwner&);
};

class PlanOwner
{
public:
    PlanOwner() : plan(NULL) {}
    ~PlanOwner() { leo2_decode_plan_destroy(plan); }
    leo2_decode_plan* plan;

private:
    PlanOwner(const PlanOwner&);
    PlanOwner& operator=(const PlanOwner&);
};

void require(bool condition, const char* message)
{
    if (!condition)
        throw std::runtime_error(message);
}

void require_result(
    leo2_result actual,
    leo2_result expected,
    const char* operation)
{
    if (actual != expected)
    {
        std::cerr << operation << ": expected " << expected
                  << ", got " << actual << std::endl;
        throw std::runtime_error(operation);
    }
}

bool all_equal(const std::vector<uint8_t>& bytes, uint8_t value)
{
    return std::find_if(bytes.begin(), bytes.end(),
        [value](uint8_t byte) { return byte != value; }) == bytes.end();
}

void fill_pattern(std::vector<uint8_t>& bytes, uint32_t seed)
{
    for (size_t i = 0; i < bytes.size(); ++i)
    {
        seed = seed * 1664525u + 1013904223u;
        bytes[i] = static_cast<uint8_t>(seed >> 24);
    }
}

void create_k1_codec_and_plan(
    leo2_context* context,
    leo2_field field,
    CodecOwner& codec,
    PlanOwner& plan)
{
    require_result(leo2_codec_create(context, 1, 1,
        LEO2_PROFILE_LEGACY_HIGH_V1, field, NULL, &codec.codec),
        LEO2_SUCCESS, "K1 codec create");
    const uint8_t original_present[1] = { 0 };
    const uint8_t recovery_present[1] = { 1 };
    require_result(leo2_decode_plan_create(codec.codec, original_present,
        recovery_present, &plan.plan), LEO2_SUCCESS, "K1 plan create");
    require(leo2_decode_plan_missing_original_count(plan.plan) == 1,
        "K1 plan lost its missing original");
}

void check_guarded_copy(
    const std::vector<uint8_t>& output,
    size_t offset,
    const std::vector<uint8_t>& expected,
    uint8_t guard,
    const char* operation)
{
    require(output.size() == expected.size() + 2 * offset,
        "guarded output has the wrong size");
    require(std::equal(expected.begin(), expected.end(),
            output.begin() + static_cast<ptrdiff_t>(offset)), operation);
    for (size_t i = 0; i < offset; ++i)
        require(output[i] == guard, "K1 copy changed a prefix guard");
    for (size_t i = offset + expected.size(); i < output.size(); ++i)
        require(output[i] == guard, "K1 copy changed a suffix guard");
}

void test_valid_single(
    leo2_context* context,
    leo2_field field,
    size_t bytes)
{
    CodecOwner codec;
    PlanOwner plan;
    create_k1_codec_and_plan(context, field, codec, plan);

    size_t scratch_bytes = 1;
    require_result(leo2_decode_plan_scratch_size(plan.plan, bytes,
        &scratch_bytes), LEO2_SUCCESS, "K1 plan scratch query");
    require(scratch_bytes == 0, "K1 plan unexpectedly requires scratch");
    scratch_bytes = 1;
    require_result(leo2_decode_scratch_size(codec.codec, bytes,
        &scratch_bytes), LEO2_SUCCESS, "K1 one-shot scratch query");
    require(scratch_bytes == 0,
        "K1 one-shot unexpectedly requires scratch");

    std::vector<uint8_t> parity(bytes);
    fill_pattern(parity, static_cast<uint32_t>(bytes * 17u + field));
    const size_t guard_bytes = 5;
    std::vector<uint8_t> output(bytes + 2 * guard_bytes, 0xa5);
    const void* original[1] = { NULL };
    const void* recovery[1] = { parity.data() };
    void* restored[1] = { output.data() + guard_bytes };

    require_result(leo2_decode_plan_execute(plan.plan, bytes, original,
        recovery, restored, NULL, 0), LEO2_SUCCESS,
        "K1 plan execute");
    check_guarded_copy(output, guard_bytes, parity, 0xa5,
        "K1 plan restored wrong bytes");

    std::fill(output.begin(), output.end(), 0x6b);
    leo2_decode_batch_item item = {
        bytes, original, recovery, restored, NULL, 0
    };
    require_result(leo2_decode_plan_execute_batch(plan.plan, &item, 1),
        LEO2_SUCCESS, "K1 one-item batch execute");
    check_guarded_copy(output, guard_bytes, parity, 0x6b,
        "K1 one-item batch restored wrong bytes");

    std::fill(output.begin(), output.end(), 0x3c);
    const uint8_t original_present[1] = { 0 };
    const uint8_t recovery_present[1] = { 1 };
    require_result(leo2_decode(codec.codec, bytes, original_present,
        recovery_present, original, recovery, restored, NULL, 0),
        LEO2_SUCCESS, "K1 one-shot execute");
    check_guarded_copy(output, guard_bytes, parity, 0x3c,
        "K1 one-shot restored wrong bytes");

    AlignedBuffer optional_scratch(64);
    std::fill(output.begin(), output.end(), 0x7d);
    require_result(leo2_decode_plan_execute(plan.plan, bytes, original,
        recovery, restored, optional_scratch.data(), optional_scratch.size()),
        LEO2_SUCCESS, "K1 optional-scratch execute");
    check_guarded_copy(output, guard_bytes, parity, 0x7d,
        "K1 optional-scratch restored wrong bytes");
}

void test_invalid_range_alias_and_atomicity(
    leo2_context* context,
    leo2_field field)
{
    const size_t bytes = field == LEO2_FIELD_GF8 ? 17 : 18;
    CodecOwner codec;
    PlanOwner plan;
    create_k1_codec_and_plan(context, field, codec, plan);

    std::vector<uint8_t> parity(bytes);
    fill_pattern(parity, 0x91e10da5u + field);
    std::vector<uint8_t> output(bytes, 0xa5);
    const void* original[1] = { NULL };
    const void* recovery[1] = { parity.data() };
    void* restored[1] = { output.data() };

    const auto require_atomic_failure = [&](leo2_result expected,
            const char* operation, const void* const* original_arg,
            const void* const* recovery_arg, void* const* restored_arg,
            void* scratch, size_t supplied_scratch) {
        std::fill(output.begin(), output.end(), 0xa5);
        require_result(leo2_decode_plan_execute(plan.plan, bytes,
            original_arg, recovery_arg, restored_arg, scratch,
            supplied_scratch), expected, operation);
        require(all_equal(output, 0xa5),
            "rejected K1 plan call modified its output");
    };

    require_atomic_failure(LEO2_INVALID_ARGUMENT, "K1 null original array",
        NULL, recovery, restored, NULL, 0);
    require_atomic_failure(LEO2_INVALID_ARGUMENT, "K1 null recovery array",
        original, NULL, restored, NULL, 0);
    require_atomic_failure(LEO2_INVALID_ARGUMENT, "K1 null restored array",
        original, recovery, NULL, NULL, 0);

    const void* unexpected_original[1] = { parity.data() };
    require_atomic_failure(LEO2_INVALID_ARGUMENT,
        "K1 nonnull missing original", unexpected_original, recovery,
        restored, NULL, 0);
    const void* null_recovery[1] = { NULL };
    require_atomic_failure(LEO2_INVALID_ARGUMENT, "K1 null recovery shard",
        original, null_recovery, restored, NULL, 0);
    void* null_restored[1] = { NULL };
    require_atomic_failure(LEO2_INVALID_ARGUMENT, "K1 null output shard",
        original, recovery, null_restored, NULL, 0);

    std::fill(output.begin(), output.end(), 0xa5);
    require_result(leo2_decode_plan_execute(plan.plan, 0, original, recovery,
        restored, NULL, 0), LEO2_INVALID_ARGUMENT, "K1 zero shard bytes");
    require(all_equal(output, 0xa5),
        "zero-byte K1 rejection modified output");

    const uintptr_t metadata_end = std::numeric_limits<uintptr_t>::max() -
        sizeof(void*) + 1;
    const void* const* const overflowing_const_metadata =
        reinterpret_cast<const void* const*>(metadata_end);
    void* const* const overflowing_mutable_metadata =
        reinterpret_cast<void* const*>(metadata_end);
    require_atomic_failure(LEO2_INVALID_ARGUMENT,
        "K1 overflowing original metadata", overflowing_const_metadata,
        recovery, restored, NULL, 0);
    require_atomic_failure(LEO2_INVALID_ARGUMENT,
        "K1 overflowing recovery metadata", original,
        overflowing_const_metadata, restored, NULL, 0);
    require_atomic_failure(LEO2_INVALID_ARGUMENT,
        "K1 overflowing restored metadata", original, recovery,
        overflowing_mutable_metadata, NULL, 0);

    const uintptr_t shard_end = std::numeric_limits<uintptr_t>::max() -
        bytes + 1;
    const void* overflowing_recovery[1] = {
        reinterpret_cast<const void*>(shard_end)
    };
    require_atomic_failure(LEO2_INVALID_ARGUMENT,
        "K1 overflowing recovery shard", original, overflowing_recovery,
        restored, NULL, 0);
    void* overflowing_output[1] = { reinterpret_cast<void*>(shard_end) };
    require_atomic_failure(LEO2_INVALID_ARGUMENT,
        "K1 overflowing output shard", original, recovery,
        overflowing_output, NULL, 0);

    void* recovery_alias_output[1] = {
        const_cast<uint8_t*>(parity.data())
    };
    require_atomic_failure(LEO2_OVERLAP, "K1 recovery/output overlap",
        original, recovery, recovery_alias_output, NULL, 0);
    void* metadata_alias_output[1];
    metadata_alias_output[0] = metadata_alias_output;
    require_atomic_failure(LEO2_OVERLAP, "K1 output/metadata overlap",
        original, recovery, metadata_alias_output, NULL, 0);

    AlignedBuffer optional_scratch(64);
    require_atomic_failure(LEO2_BAD_ALIGNMENT,
        "K1 misaligned optional scratch", original, recovery, restored,
        static_cast<uint8_t*>(optional_scratch.data()) + 1, 63);
    void* scratch_alias_output[1] = { optional_scratch.data() };
    const std::vector<uint8_t> scratch_before(64, 0);
    require_atomic_failure(LEO2_OVERLAP, "K1 output/scratch overlap",
        original, recovery, scratch_alias_output, optional_scratch.data(),
        optional_scratch.size());
    require(memcmp(optional_scratch.data(), scratch_before.data(), 64) == 0,
        "rejected K1 scratch overlap modified scratch");

    const uint8_t missing_original[1] = { 0 };
    const uint8_t received_recovery[1] = { 1 };
    uint8_t invalid_original[1] = { 2 };
    uint8_t invalid_recovery[1] = { 2 };
    std::fill(output.begin(), output.end(), 0xa5);
    require_result(leo2_decode(codec.codec, bytes, invalid_original,
        received_recovery, original, recovery, restored, NULL, 0),
        LEO2_INVALID_ARGUMENT, "K1 invalid original presence");
    require_result(leo2_decode(codec.codec, bytes, missing_original,
        invalid_recovery, original, recovery, restored, NULL, 0),
        LEO2_INVALID_ARGUMENT, "K1 invalid recovery presence");
    require_result(leo2_decode(codec.codec, bytes, NULL,
        received_recovery, original, recovery, restored, NULL, 0),
        LEO2_INVALID_ARGUMENT, "K1 null original presence");
    require_result(leo2_decode(codec.codec, bytes, missing_original,
        NULL, original, recovery, restored, NULL, 0),
        LEO2_INVALID_ARGUMENT, "K1 null recovery presence");
    require_result(leo2_decode(codec.codec, bytes,
        reinterpret_cast<const uint8_t*>(
            std::numeric_limits<uintptr_t>::max()),
        received_recovery, original, recovery, restored, NULL, 0),
        LEO2_INVALID_ARGUMENT, "K1 overflowing original presence");
    require_result(leo2_decode(codec.codec, bytes, missing_original,
        reinterpret_cast<const uint8_t*>(
            std::numeric_limits<uintptr_t>::max()),
        original, recovery, restored, NULL, 0),
        LEO2_INVALID_ARGUMENT, "K1 overflowing recovery presence");
    require(all_equal(output, 0xa5),
        "invalid K1 one-shot presence modified output");

    std::vector<uint8_t> original_presence(bytes, 0);
    void* original_presence_output[1] = { original_presence.data() };
    const std::vector<uint8_t> original_presence_before = original_presence;
    require_result(leo2_decode(codec.codec, bytes,
        original_presence.data(), received_recovery, original, recovery,
        original_presence_output, NULL, 0), LEO2_OVERLAP,
        "K1 output/original-presence overlap");
    require(original_presence == original_presence_before,
        "K1 original-presence overlap modified metadata");

    std::vector<uint8_t> recovery_presence(bytes, 0);
    recovery_presence[0] = 1;
    void* recovery_presence_output[1] = { recovery_presence.data() };
    const std::vector<uint8_t> recovery_presence_before = recovery_presence;
    require_result(leo2_decode(codec.codec, bytes, missing_original,
        recovery_presence.data(), original, recovery,
        recovery_presence_output, NULL, 0), LEO2_OVERLAP,
        "K1 output/recovery-presence overlap");
    require(recovery_presence == recovery_presence_before,
        "K1 recovery-presence overlap modified metadata");

    alignas(64) uint8_t presence_scratch[64];
    memset(presence_scratch, 0, sizeof(presence_scratch));
    const std::vector<uint8_t> output_before = output;
    require_result(leo2_decode(codec.codec, bytes, presence_scratch,
        received_recovery, original, recovery, restored, presence_scratch,
        sizeof(presence_scratch)), LEO2_OVERLAP,
        "K1 scratch/original-presence overlap");
    require(output == output_before,
        "K1 presence/scratch rejection modified output");
    for (size_t i = 0; i < sizeof(presence_scratch); ++i)
        require(presence_scratch[i] == 0,
            "K1 presence/scratch rejection modified metadata");

    if (field == LEO2_FIELD_GF16)
    {
        std::fill(output.begin(), output.end(), 0xa5);
        require_result(leo2_decode_plan_execute(plan.plan, bytes - 1,
            original, recovery, restored, NULL, 0), LEO2_UNSUPPORTED,
            "K1 GF16 odd physical length");
        require_result(leo2_decode(codec.codec, bytes - 1, missing_original,
            received_recovery, original, recovery, restored, NULL, 0),
            LEO2_UNSUPPORTED, "K1 GF16 one-shot odd physical length");
        require(all_equal(output, 0xa5),
            "unsupported K1 GF16 tail modified output");
    }
}

void test_batch_paths(
    leo2_context* context,
    leo2_field field)
{
    const size_t bytes = field == LEO2_FIELD_GF8 ? 65 : 66;
    const size_t item_count = 9;
    CodecOwner codec;
    PlanOwner plan;
    create_k1_codec_and_plan(context, field, codec, plan);

    std::vector<std::vector<uint8_t> > parity(
        item_count, std::vector<uint8_t>(bytes));
    std::vector<std::vector<uint8_t> > output(
        item_count, std::vector<uint8_t>(bytes, 0xa5));
    std::vector<std::array<const void*, 1> > original(item_count);
    std::vector<std::array<const void*, 1> > recovery(item_count);
    std::vector<std::array<void*, 1> > restored(item_count);
    std::vector<leo2_decode_batch_item> items(item_count);
    for (size_t i = 0; i < item_count; ++i)
    {
        fill_pattern(parity[i], static_cast<uint32_t>(0x9e3779b9u + i));
        original[i][0] = NULL;
        recovery[i][0] = parity[i].data();
        restored[i][0] = output[i].data();
        items[i].shard_bytes = bytes;
        items[i].original = original[i].data();
        items[i].recovery = recovery[i].data();
        items[i].restored_original = restored[i].data();
        items[i].scratch = NULL;
        items[i].scratch_bytes = 0;
    }

    require_result(leo2_decode_plan_execute_batch(plan.plan, items.data(), 2),
        LEO2_SUCCESS, "K1 ordinary multi-item batch");
    require(output[0] == parity[0] && output[1] == parity[1],
        "K1 ordinary batch restored wrong bytes");

    size_t preflight_bytes = 0;
    require_result(leo2_decode_plan_batch_preflight_scratch_size(
        plan.plan, item_count, &preflight_bytes), LEO2_SUCCESS,
        "K1 scalable preflight scratch query");
    require(preflight_bytes != 0 &&
            preflight_bytes % leo2_scratch_alignment() == 0,
        "K1 scalable preflight scratch is missing or unaligned");
    AlignedBuffer preflight(preflight_bytes);

    for (size_t i = 0; i < item_count; ++i)
        std::fill(output[i].begin(), output[i].end(), 0x6b);
    require_result(leo2_decode_plan_execute_batch_with_preflight_scratch(
        plan.plan, items.data(), items.size(), preflight.data(),
        preflight.size()), LEO2_SUCCESS, "K1 scalable batch execute");
    for (size_t i = 0; i < item_count; ++i)
        require(output[i] == parity[i],
            "K1 scalable batch restored wrong bytes");

    for (size_t i = 0; i < item_count; ++i)
        std::fill(output[i].begin(), output[i].end(), 0x7d);
    const std::vector<std::vector<uint8_t> > failure_before = output;
    const void* const* const saved_recovery = items.back().recovery;
    items.back().recovery = NULL;
    require_result(leo2_decode_plan_execute_batch_with_preflight_scratch(
        plan.plan, items.data(), items.size(), preflight.data(),
        preflight.size()), LEO2_INVALID_ARGUMENT,
        "K1 scalable late-item invalid metadata");
    require(output == failure_before,
        "K1 scalable validation failure partially restored outputs");
    items.back().recovery = saved_recovery;

    const std::vector<uint8_t> first_parity_before = parity[0];
    void* const saved_output = restored.back()[0];
    restored.back()[0] = parity[0].data();
    require_result(leo2_decode_plan_execute_batch_with_preflight_scratch(
        plan.plan, items.data(), items.size(), preflight.data(),
        preflight.size()), LEO2_OVERLAP,
        "K1 scalable cross-item output/input overlap");
    require(output == failure_before && parity[0] == first_parity_before,
        "K1 scalable alias rejection modified caller storage");
    restored.back()[0] = saved_output;

    void* const* const saved_restored_metadata =
        items.back().restored_original;
    items.back().restored_original = reinterpret_cast<void* const*>(
        std::numeric_limits<uintptr_t>::max() - sizeof(void*) + 1);
    require_result(leo2_decode_plan_execute_batch_with_preflight_scratch(
        plan.plan, items.data(), items.size(), preflight.data(),
        preflight.size()), LEO2_INVALID_ARGUMENT,
        "K1 scalable overflowing late-item metadata");
    require(output == failure_before,
        "K1 scalable range rejection partially restored outputs");
    items.back().restored_original = saved_restored_metadata;
}

void test_field(
    leo2_context* context,
    leo2_field field)
{
    if (field == LEO2_FIELD_GF8)
    {
        static const size_t sizes[] = { 1, 17, 65, 257 };
        for (size_t i = 0; i < sizeof(sizes) / sizeof(sizes[0]); ++i)
            test_valid_single(context, field, sizes[i]);
    }
    else
    {
        static const size_t sizes[] = { 2, 18, 66, 258 };
        for (size_t i = 0; i < sizeof(sizes) / sizeof(sizes[0]); ++i)
            test_valid_single(context, field, sizes[i]);
    }
    test_invalid_range_alias_and_atomicity(context, field);
    test_batch_paths(context, field);
}

} // namespace

int main()
{
    try
    {
#if defined(LEO2_EXPECT_K1_DECODE_VALIDATOR_INLINE)
        const bool expected_inline = true;
#else
        const bool expected_inline = false;
#endif
        require(leopard2_internal::
                K1DecodeValidatorInlineExperimentEnabled() == expected_inline,
            "K1 decode-validator inline experiment marker mismatch");

        leo2_context_options options;
        memset(&options, 0, sizeof(options));
        options.struct_size = sizeof(options);
        options.backend = LEO2_BACKEND_AUTO;
        options.thread_count = 1;
        leo2_context* context = NULL;
        require_result(leo2_context_create(&options, &context),
            LEO2_SUCCESS, "K1 production context create");
        const uint32_t field_mask = leo2_context_field_mask(context);
        require(field_mask != 0, "K1 production context has no field");
        if ((field_mask & LEO2_FIELD_MASK_GF8) != 0)
            test_field(context, LEO2_FIELD_GF8);
        if ((field_mask & LEO2_FIELD_MASK_GF16) != 0)
            test_field(context, LEO2_FIELD_GF16);
        leo2_context_destroy(context);

        std::cout << "{\"k1_inline\":" << (expected_inline ? 1 : 0)
                  << ",\"production_behavior\":\"ok\"}" << std::endl;
        return 0;
    }
    catch (const std::exception& error)
    {
        std::cerr << error.what() << std::endl;
        return 1;
    }
}
