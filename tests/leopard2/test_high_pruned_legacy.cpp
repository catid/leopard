/*
    Copyright (c) 2017 Christopher A. Taylor.  All rights reserved.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the conditions in the Leopard-RS
    BSD license are met.  See the repository LICENSE file.
*/

/*
    Exact Leopard1 diagnostic for the production Algorithm 5 schedule path.
    Old and new encoders must emit identical parity; old Leopard's full decoder
    and both new pruned workspaces must restore identical original bytes.
*/

#include "leopard.h"
#include "leopard2.h"
#if defined(LEO2_REQUIRE_HIGH_PRUNED_STAGE_HOOKS) && \
    !defined(LEO2_ENABLE_TEST_HOOKS)
#error "high-pruned stage counter target requires Leopard2 test hooks"
#endif
#if defined(LEO2_ENABLE_TEST_HOOKS)
#include "Leopard2Backend.h"
#include "Leopard2Plan.h"
#include "LeopardFF16.h"
#include "LeopardFF8.h"
#endif

#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <new>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#if defined(_MSC_VER)
#include <malloc.h>
#endif

namespace {

class AlignedBytes
{
public:
    explicit AlignedBytes(size_t bytes)
        : data_(NULL), bytes_(bytes)
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

    ~AlignedBytes()
    {
#if defined(_MSC_VER)
        _aligned_free(data_);
#else
        free(data_);
#endif
    }

    uint8_t* data() { return static_cast<uint8_t*>(data_); }
    const uint8_t* data() const { return static_cast<const uint8_t*>(data_); }
    size_t size() const { return bytes_; }

private:
    AlignedBytes(const AlignedBytes&);
    AlignedBytes& operator=(const AlignedBytes&);
    void* data_;
    size_t bytes_;
};

void require(bool condition, const char* message)
{
    if (!condition)
        throw std::runtime_error(message);
}

void require_result(leo2_result result, const char* operation)
{
    if (result != LEO2_SUCCESS)
        throw std::runtime_error(std::string(operation) + ": " +
            leo2_result_string(result));
}

uint32_t mix(uint32_t value)
{
    value ^= value >> 16;
    value *= 0x7feb352du;
    value ^= value >> 15;
    value *= 0x846ca68bu;
    return value ^ (value >> 16);
}

struct Case
{
    uint32_t k;
    uint32_t r;
    leo2_field field;
    size_t bytes;
};

#if defined(LEO2_ENABLE_TEST_HOOKS)
void verify_gf8_one_loss_receive_stage(
    leo2_context* context,
    const leo2_codec* codec,
    const Case& test,
    const std::vector<const void*>& original_input,
    const std::vector<void*>& recovery_output)
{
    if (test.field != LEO2_FIELD_GF8 || test.k != 240 || test.r != 16 ||
        test.bytes != 64)
        return;

    std::vector<uint8_t> original_present(test.k, 1);
    std::vector<uint8_t> recovery_present(test.r, 1);
    original_present[0] = 0;
    leo2_decode_plan* plan = NULL;
    require_result(leo2_decode_plan_create(codec, &original_present[0],
        &recovery_present[0], &plan), "one-loss plan create");

    size_t scratch_bytes = 0;
    require_result(leo2_decode_plan_scratch_size(
        plan, test.bytes, &scratch_bytes), "one-loss scratch query");
    AlignedBytes scratch(scratch_bytes);
    AlignedBytes restored_storage((test.bytes + 63u) & ~size_t(63u));
    std::vector<const void*> original = original_input;
    original[0] = NULL;
    std::vector<const void*> recovery(test.r);
    for (uint32_t i = 0; i < test.r; ++i)
        recovery[i] = recovery_output[i];
    std::vector<void*> restored(test.k, NULL);
    restored[0] = restored_storage.data();

    leopard::ff8::TestOnlyResetHighDecodeCounts();
    require_result(leo2_decode_plan_execute(plan, test.bytes,
        &original[0], &recovery[0], &restored[0],
        scratch.data(), scratch.size()), "one-loss decode");
    const leopard::ff8::TestOnlyHighDecodeCounts counts =
        leopard::ff8::TestOnlyGetHighDecodeCounts();
    std::cout << "COUNTERS high_receive_gf8_k240_r16_l1 backend="
              << static_cast<unsigned>(leo2_context_backend(context))
              << " fused4=" << counts.receive_ifft_butterfly4_out_of_place
              << " copies=" << counts.receive_copy_shards
              << " zeros=" << counts.receive_zero_shards
              << " block0_ifft_elided="
              << counts.syndrome_block_zero_ifft_elisions
              << " syndrome_fft=" << counts.syndrome_forward_transforms
              << " block0_xor_shards="
              << counts.syndrome_block_zero_xor_shards << std::endl;
    require(memcmp(restored[0], original_input[0], test.bytes) == 0,
        "one-loss restored original mismatch");

    const leo2_backend selected = leo2_context_backend(context);
    if (selected == LEO2_BACKEND_SSSE3 ||
        selected == LEO2_BACKEND_AVX2 ||
        selected == LEO2_BACKEND_AVX512)
    {
        require(counts.receive_ifft_butterfly4_out_of_place == 56,
            "qualified GF8 one-loss fused-group count drifted");
        require(counts.receive_copy_shards == 15,
            "qualified GF8 one-loss staged-copy count drifted");
        require(counts.receive_zero_shards == 1,
            "qualified GF8 one-loss staged-zero count drifted");
    }
    else
    {
        require(counts.receive_ifft_butterfly4_out_of_place == 0,
            "unqualified GF8 backend fused receive rows");
        require(counts.receive_copy_shards == test.k - 1,
            "copy-first GF8 one-loss staged-copy count drifted");
        require(counts.receive_zero_shards == 1,
            "copy-first GF8 one-loss staged-zero count drifted");
    }
    require(counts.syndrome_block_zero_ifft_elisions == 1,
        "GF8 block-zero inverse transform was not elided");
    require(counts.syndrome_forward_transforms == 1 &&
            counts.syndrome_forward_transform_elisions == 0,
        "GF8 later-block syndrome did not retain exactly one forward transform");
    require(counts.syndrome_block_zero_xor_shards == 1,
        "GF8 block-zero raw-source XOR count drifted");
    leo2_decode_plan_destroy(plan);
}

template<class Counts, class ResetCounts, class GetCounts>
void verify_block_zero_only_cancellation_case(
    leo2_context* context,
    leo2_field field,
    uint32_t k,
    size_t bytes,
    uint32_t layout_flag,
    ResetCounts reset_counts,
    GetCounts get_counts,
    const char* label)
{
    leo2_codec_options specialized_options = {};
    specialized_options.struct_size = sizeof(specialized_options);
    specialized_options.flags = LEO2_CODEC_FORCE_SPECIALIZED_DECODE |
        layout_flag;
    leo2_codec* specialized_codec = NULL;
    require_result(leo2_codec_create(context, k, k,
        LEO2_PROFILE_LEGACY_HIGH_V1, field, &specialized_options,
        &specialized_codec), "block-zero-only specialized codec create");

    leo2_codec_options generic_options = {};
    generic_options.struct_size = sizeof(generic_options);
    generic_options.flags = LEO2_CODEC_FORCE_GENERIC_DECODE;
    leo2_codec* generic_codec = NULL;
    require_result(leo2_codec_create(context, k, k,
        LEO2_PROFILE_LEGACY_HIGH_V1, field, &generic_options,
        &generic_codec), "block-zero-only generic codec create");

    const size_t stride = (bytes + 63u) & ~size_t(63u);
    AlignedBytes original_storage(static_cast<size_t>(k) * stride);
    AlignedBytes recovery_storage(static_cast<size_t>(k) * stride);
    std::vector<const void*> originals(k);
    std::vector<void*> recovery_output(k);
    for (uint32_t shard = 0; shard < k; ++shard)
    {
        uint8_t* original = original_storage.data() +
            static_cast<size_t>(shard) * stride;
        for (size_t byte = 0; byte < bytes; ++byte)
            original[byte] = static_cast<uint8_t>(
                mix(shard * 65537u + static_cast<uint32_t>(byte) * 17u));
        originals[shard] = original;
        recovery_output[shard] = recovery_storage.data() +
            static_cast<size_t>(shard) * stride;
    }
    size_t encode_scratch_bytes = 0;
    require_result(leo2_encode_scratch_size(
        specialized_codec, bytes, &encode_scratch_bytes),
        "block-zero-only encode scratch query");
    AlignedBytes encode_scratch(encode_scratch_bytes);
    require_result(leo2_encode(specialized_codec, bytes, &originals[0],
        &recovery_output[0], encode_scratch.data(), encode_scratch.size()),
        "block-zero-only encode");

    std::vector<uint8_t> original_present(k, 0);
    std::vector<uint8_t> recovery_present(k, 1);
    leo2_decode_plan* specialized_plan = NULL;
    leo2_decode_plan* generic_plan = NULL;
    require_result(leo2_decode_plan_create(specialized_codec,
        &original_present[0], &recovery_present[0], &specialized_plan),
        "block-zero-only specialized plan create");
    require_result(leo2_decode_plan_create(generic_codec,
        &original_present[0], &recovery_present[0], &generic_plan),
        "block-zero-only generic plan create");

    std::vector<const void*> missing_originals(k, NULL);
    std::vector<const void*> recovery(k);
    for (uint32_t i = 0; i < k; ++i)
        recovery[i] = recovery_output[i];
    AlignedBytes specialized_storage(static_cast<size_t>(k) * stride);
    AlignedBytes generic_storage(static_cast<size_t>(k) * stride);
    std::vector<void*> specialized_output(k), generic_output(k);
    for (uint32_t i = 0; i < k; ++i)
    {
        specialized_output[i] = specialized_storage.data() +
            static_cast<size_t>(i) * stride;
        generic_output[i] = generic_storage.data() +
            static_cast<size_t>(i) * stride;
    }

    size_t specialized_scratch_bytes = 0;
    size_t generic_scratch_bytes = 0;
    require_result(leo2_decode_plan_scratch_size(
        specialized_plan, bytes, &specialized_scratch_bytes),
        "block-zero-only specialized scratch query");
    require_result(leo2_decode_plan_scratch_size(
        generic_plan, bytes, &generic_scratch_bytes),
        "block-zero-only generic scratch query");
    AlignedBytes specialized_scratch(specialized_scratch_bytes);
    AlignedBytes generic_scratch(generic_scratch_bytes);

    reset_counts();
    require_result(leo2_decode_plan_execute(specialized_plan, bytes,
        &missing_originals[0], &recovery[0], &specialized_output[0],
        specialized_scratch.data(), specialized_scratch.size()),
        "block-zero-only specialized decode");
    const Counts counts = get_counts();
    require_result(leo2_decode_plan_execute(generic_plan, bytes,
        &missing_originals[0], &recovery[0], &generic_output[0],
        generic_scratch.data(), generic_scratch.size()),
        "block-zero-only generic decode");
    for (uint32_t i = 0; i < k; ++i)
    {
        require(memcmp(specialized_output[i], originals[i], bytes) == 0,
            "block-zero-only specialized recovery mismatch");
        require(memcmp(generic_output[i], specialized_output[i], bytes) == 0,
            "block-zero-only specialized/generic mismatch");
    }

    const uint64_t kernel_passes = bytes > 64 && (bytes & 63u) != 0 ? 2 : 1;
    require(counts.syndrome_block_zero_ifft_elisions == kernel_passes,
        "block-zero-only inverse-elision count mismatch");
    require(counts.syndrome_forward_transforms == 0,
        "block-zero-only path retained a forward syndrome transform");
    require(counts.syndrome_forward_transform_elisions == kernel_passes,
        "block-zero-only forward-elision count mismatch");
    require(counts.syndrome_block_zero_xor_shards == 0,
        "block-zero-only path unexpectedly XORed staged raw sources");
    std::cout << "COUNTERS " << label
              << " block0_ifft_elided="
              << counts.syndrome_block_zero_ifft_elisions
              << " syndrome_fft=" << counts.syndrome_forward_transforms
              << " syndrome_fft_elided="
              << counts.syndrome_forward_transform_elisions << std::endl;

    leo2_decode_plan_destroy(generic_plan);
    leo2_decode_plan_destroy(specialized_plan);
    leo2_codec_destroy(generic_codec);
    leo2_codec_destroy(specialized_codec);
}

void verify_block_zero_only_cancellation(leo2_context* context)
{
    verify_block_zero_only_cancellation_case<
        leopard::ff8::TestOnlyHighDecodeCounts>(
        context, LEO2_FIELD_GF8, 128, 65,
        LEO2_CODEC_FORCE_MATERIALIZED_DECODE,
        leopard::ff8::TestOnlyResetHighDecodeCounts,
        leopard::ff8::TestOnlyGetHighDecodeCounts,
        "alg5_block0_only_gf8_materialized");
    verify_block_zero_only_cancellation_case<
        leopard::ff8::TestOnlyHighDecodeCounts>(
        context, LEO2_FIELD_GF8, 128, 64,
        LEO2_CODEC_FORCE_TILED_DECODE,
        leopard::ff8::TestOnlyResetHighDecodeCounts,
        leopard::ff8::TestOnlyGetHighDecodeCounts,
        "alg5_block0_only_gf8_tiled");
    verify_block_zero_only_cancellation_case<
        leopard::ff16::TestOnlyHighDecodeCounts>(
        context, LEO2_FIELD_GF16, 256, 66,
        LEO2_CODEC_FORCE_MATERIALIZED_DECODE,
        leopard::ff16::TestOnlyResetHighDecodeCounts,
        leopard::ff16::TestOnlyGetHighDecodeCounts,
        "alg5_block0_only_gf16_materialized");
    verify_block_zero_only_cancellation_case<
        leopard::ff16::TestOnlyHighDecodeCounts>(
        context, LEO2_FIELD_GF16, 256, 64,
        LEO2_CODEC_FORCE_TILED_DECODE,
        leopard::ff16::TestOnlyResetHighDecodeCounts,
        leopard::ff16::TestOnlyGetHighDecodeCounts,
        "alg5_block0_only_gf16_tiled");
}

struct GF8CancellationField
{
    typedef leopard::ff8::ffe_t Ffe;
    typedef leopard::ff8::TestOnlyHighDecodeCounts Counts;
    static const unsigned kModulus = leopard::ff8::kModulus;

    static void prepare(unsigned n, unsigned side, Ffe* factors)
    {
        leopard::ff8::PrepareHighDecode(n, side, factors);
    }
    static bool prepare_pruned(unsigned size, unsigned shift, bool inverse,
        const uint8_t* input_mask, const uint8_t* output_mask,
        leopard2_internal::PrunedTransformPlan& plan)
    {
        return leopard::ff8::PreparePrunedTransformPlan(
            size, shift, inverse, input_mask, output_mask, plan);
    }
    static void prepared(uint64_t bytes, unsigned n, unsigned side,
        const void* const* inputs, const uint8_t* requested,
        const Ffe* locator, const Ffe* factors, void** work)
    {
        leopard::ff8::ReedSolomonDecodeHighPrepared(
            bytes, n, side, inputs, requested, locator, factors, work);
    }
    static void pruned(uint64_t bytes, unsigned n, unsigned side,
        const void* const* inputs, const uint16_t* prefixes,
        const uint32_t* requested,
        const leopard2_internal::DecodeOutputBlock* output_blocks,
        unsigned output_block_count, const Ffe* locator, const Ffe* factors,
        const leopard2_internal::PrunedTransformBlock* input_plans,
        unsigned input_plan_count, void** work)
    {
        leopard::ff8::ReedSolomonDecodeHighPrunedPlanned(
            leopard::backend::GetDefaultOps(), bytes, n, side, inputs,
            prefixes, requested, output_blocks, output_block_count, locator,
            factors, input_plans, input_plan_count, NULL, 0, NULL, work);
    }
    static void tiled_pruned(uint64_t bytes, unsigned n, unsigned side,
        const void* const* inputs, const uint16_t* prefixes,
        const uint32_t* requested,
        const leopard2_internal::DecodeOutputBlock* output_blocks,
        unsigned output_block_count, const Ffe* locator, const Ffe* factors,
        void* const* outputs,
        const leopard2_internal::PrunedTransformBlock* input_plans,
        unsigned input_plan_count, void** work)
    {
        leopard::ff8::ReedSolomonDecodeHighTiledPrunedPlanned(
            leopard::backend::GetDefaultOps(), bytes, n, side, inputs,
            prefixes, requested, output_blocks, output_block_count, locator,
            factors, outputs, input_plans, input_plan_count, NULL, 0, work);
    }
    static void reset() { leopard::ff8::TestOnlyResetHighDecodeCounts(); }
    static Counts counts() { return leopard::ff8::TestOnlyGetHighDecodeCounts(); }
};

struct GF16CancellationField
{
    typedef leopard::ff16::ffe_t Ffe;
    typedef leopard::ff16::TestOnlyHighDecodeCounts Counts;
    static const unsigned kModulus = leopard::ff16::kModulus;

    static void prepare(unsigned n, unsigned side, Ffe* factors)
    {
        leopard::ff16::PrepareHighDecode(n, side, factors);
    }
    static bool prepare_pruned(unsigned size, unsigned shift, bool inverse,
        const uint8_t* input_mask, const uint8_t* output_mask,
        leopard2_internal::PrunedTransformPlan& plan)
    {
        return leopard::ff16::PreparePrunedTransformPlan(
            size, shift, inverse, input_mask, output_mask, plan);
    }
    static void prepared(uint64_t bytes, unsigned n, unsigned side,
        const void* const* inputs, const uint8_t* requested,
        const Ffe* locator, const Ffe* factors, void** work)
    {
        leopard::ff16::ReedSolomonDecodeHighPrepared(
            bytes, n, side, inputs, requested, locator, factors, work);
    }
    static void pruned(uint64_t bytes, unsigned n, unsigned side,
        const void* const* inputs, const uint16_t* prefixes,
        const uint32_t* requested,
        const leopard2_internal::DecodeOutputBlock* output_blocks,
        unsigned output_block_count, const Ffe* locator, const Ffe* factors,
        const leopard2_internal::PrunedTransformBlock* input_plans,
        unsigned input_plan_count, void** work)
    {
        leopard::ff16::ReedSolomonDecodeHighPrunedPlanned(
            leopard::backend::GetDefaultOps(), bytes, n, side, inputs,
            prefixes, requested, output_blocks, output_block_count, locator,
            factors, input_plans, input_plan_count, NULL, 0, NULL, work);
    }
    static void tiled_pruned(uint64_t bytes, unsigned n, unsigned side,
        const void* const* inputs, const uint16_t* prefixes,
        const uint32_t* requested,
        const leopard2_internal::DecodeOutputBlock* output_blocks,
        unsigned output_block_count, const Ffe* locator, const Ffe* factors,
        void* const* outputs,
        const leopard2_internal::PrunedTransformBlock* input_plans,
        unsigned input_plan_count, void** work)
    {
        leopard::ff16::ReedSolomonDecodeHighTiledPrunedPlanned(
            leopard::backend::GetDefaultOps(), bytes, n, side, inputs,
            prefixes, requested, output_blocks, output_block_count, locator,
            factors, outputs, input_plans, input_plan_count, NULL, 0, work);
    }
    static void reset() { leopard::ff16::TestOnlyResetHighDecodeCounts(); }
    static Counts counts() { return leopard::ff16::TestOnlyGetHighDecodeCounts(); }
};

template<class Field>
void verify_internal_cancellation_pattern(
    unsigned side, bool block_zero_live, const char* label)
{
    typedef typename Field::Ffe Ffe;
    typedef typename Field::Counts Counts;
    const unsigned n = side * 8;
    const size_t bytes = 64;
    AlignedBytes sources(static_cast<size_t>(n) * bytes);
    std::vector<const void*> inputs(n, NULL);
    for (unsigned coordinate = 0; coordinate < n; ++coordinate)
    {
        uint8_t* shard = sources.data() +
            static_cast<size_t>(coordinate) * bytes;
        for (size_t byte = 0; byte < bytes; ++byte)
            shard[byte] = static_cast<uint8_t>(
                mix(coordinate * 257u + static_cast<unsigned>(byte) * 13u));
    }
    unsigned block_zero_live_rows = 0;
    if (block_zero_live)
    {
        inputs[0] = sources.data();
        inputs[side - 1] = sources.data() +
            static_cast<size_t>(side - 1) * bytes;
        block_zero_live_rows = 2;
    }
    // Block one is empty.  The first later contribution is sparse and owns
    // an exact input plan; the next block is dense and must accumulate into
    // the first.  Both are subsequently reused as output destinations.
    const unsigned sparse_offset = side * 2;
    inputs[sparse_offset] = sources.data() +
        static_cast<size_t>(sparse_offset) * bytes;
    inputs[sparse_offset + side / 2] = sources.data() +
        static_cast<size_t>(sparse_offset + side / 2) * bytes;
    inputs[sparse_offset + side - 1] = sources.data() +
        static_cast<size_t>(sparse_offset + side - 1) * bytes;
    const unsigned dense_offset = side * 3;
    for (unsigned i = 0; i < side; ++i)
        inputs[dense_offset + i] = sources.data() +
            static_cast<size_t>(dense_offset + i) * bytes;

    std::vector<uint16_t> prefixes(n / side, 0);
    for (unsigned coordinate = 0; coordinate < n; ++coordinate)
        if (inputs[coordinate])
            prefixes[coordinate / side] = static_cast<uint16_t>(
                coordinate % side + 1);
    const uint32_t requested_coordinates[] = {
        sparse_offset + 1, dense_offset + side / 2
    };
    const leopard2_internal::DecodeOutputBlock output_blocks[] = {
        { 2, 2, 0, 1 },
        { 3, side / 2 + 1, 1, 2 }
    };
    std::vector<uint8_t> requested_mask(n, 0);
    requested_mask[requested_coordinates[0]] = 1;
    requested_mask[requested_coordinates[1]] = 1;
    std::vector<Ffe> locator(n);
    std::vector<Ffe> factors(n);
    for (unsigned i = 0; i < n; ++i)
        locator[i] = static_cast<Ffe>(mix(i + n) % Field::kModulus);
    Field::prepare(n, side, &factors[0]);

    std::vector<leopard2_internal::PrunedTransformBlock> plans;
    std::vector<uint8_t> input_mask(side), output_mask(side, 1);
    for (unsigned block = 0; block < n / side; ++block)
    {
        unsigned live_count = 0;
        for (unsigned i = 0; i < side; ++i)
        {
            input_mask[i] = inputs[block * side + i] != NULL;
            live_count += input_mask[i];
        }
        if (live_count == 0 || live_count == side)
            continue;
        leopard2_internal::PrunedTransformBlock entry;
        entry.block = block;
        require(Field::prepare_pruned(side, block * side, true,
            &input_mask[0], &output_mask[0], entry.plan),
            "cancellation input plan build");
        plans.push_back(std::move(entry));
    }
    require(!plans.empty() &&
            plans[block_zero_live ? 1u : 0u].block == 2,
        "first later cancellation input plan is not exact-pruned");

    AlignedBytes expected_storage(static_cast<size_t>(n) * bytes);
    AlignedBytes materialized_storage(static_cast<size_t>(n) * bytes);
    std::vector<void*> expected_work(n), materialized_work(n);
    for (unsigned i = 0; i < n; ++i)
    {
        expected_work[i] = expected_storage.data() +
            static_cast<size_t>(i) * bytes;
        materialized_work[i] = materialized_storage.data() +
            static_cast<size_t>(i) * bytes;
    }
    Field::prepared(bytes, n, side, &inputs[0], &requested_mask[0],
        &locator[0], &factors[0], &expected_work[0]);

    const auto verify_counts = [block_zero_live, block_zero_live_rows, label](
        const Counts& counts) {
        require(counts.syndrome_block_zero_ifft_elisions ==
                (block_zero_live ? 1u : 0u),
            "synthetic block-zero inverse-elision count mismatch");
        require(counts.syndrome_forward_transforms == 1 &&
                counts.syndrome_forward_transform_elisions == 0,
            "synthetic later contribution forward-transform count mismatch");
        require(counts.syndrome_block_zero_xor_shards ==
                block_zero_live_rows,
            "synthetic block-zero source XOR count mismatch");
        require(counts.syndrome_materialized_blocks != 0,
            "synthetic first exact-pruned contribution was not materialized");
        require(counts.syndrome_accumulated_blocks != 0,
            "synthetic second dense contribution was not accumulated");
        std::cout << "COUNTERS " << label
                  << " block0_ifft_elided="
                  << counts.syndrome_block_zero_ifft_elisions
                  << " syndrome_materialized="
                  << counts.syndrome_materialized_blocks
                  << " syndrome_accumulated="
                  << counts.syndrome_accumulated_blocks << std::endl;
    };

    Field::reset();
    Field::pruned(bytes, n, side, &inputs[0], &prefixes[0],
        requested_coordinates, output_blocks, 2, &locator[0], &factors[0],
        &plans[0], plans.size(), &materialized_work[0]);
    Counts counts = Field::counts();
    verify_counts(counts);
    for (unsigned i = 0; i < 2; ++i)
        require(memcmp(expected_work[requested_coordinates[i]],
                materialized_work[requested_coordinates[i]], bytes) == 0,
            "materialized cancellation overwrote a former input temp incorrectly");

    AlignedBytes tiled_storage(static_cast<size_t>(side * 2 + 2) * bytes);
    std::vector<void*> tiled_work(side * 2 + 2);
    for (unsigned i = 0; i < side * 2 + 2; ++i)
        tiled_work[i] = tiled_storage.data() +
            static_cast<size_t>(i) * bytes;
    Field::reset();
    Field::tiled_pruned(bytes, n, side, &inputs[0], &prefixes[0],
        requested_coordinates, output_blocks, 2, &locator[0], &factors[0],
        &tiled_work[side * 2], &plans[0], plans.size(), &tiled_work[0]);
    counts = Field::counts();
    verify_counts(counts);
    for (unsigned i = 0; i < 2; ++i)
        require(memcmp(expected_work[requested_coordinates[i]],
                tiled_work[side * 2 + i], bytes) == 0,
            "tiled cancellation overwrote its former input tile incorrectly");
}

void verify_internal_cancellation_patterns()
{
    verify_internal_cancellation_pattern<GF8CancellationField>(
        16, false, "alg5_gf8_block0_empty");
    verify_internal_cancellation_pattern<GF8CancellationField>(
        16, true, "alg5_gf8_block0_sparse");
    verify_internal_cancellation_pattern<GF16CancellationField>(
        32, false, "alg5_gf16_block0_empty");
    verify_internal_cancellation_pattern<GF16CancellationField>(
        32, true, "alg5_gf16_block0_sparse");
}
#endif

void run_case(leo2_context* context, const Case& test, size_t case_index,
    uint64_t& parity_bytes, uint64_t& restored_bytes)
{
    const size_t stride = (test.bytes + 63u) & ~size_t(63u);
    AlignedBytes originals(static_cast<size_t>(test.k) * stride);
    std::vector<const void*> original_input(test.k);
    for (uint32_t shard = 0; shard < test.k; ++shard)
    {
        uint8_t* data = originals.data() + static_cast<size_t>(shard) * stride;
        for (size_t byte = 0; byte < test.bytes; ++byte)
            data[byte] = static_cast<uint8_t>(
                mix(static_cast<uint32_t>(case_index * 65537u +
                    shard * 257u + byte)));
        original_input[shard] = data;
    }

    const unsigned old_encode_count = leo_encode_work_count(test.k, test.r);
    require(old_encode_count != 0, "legacy encode work count");
    AlignedBytes old_encode_storage(
        static_cast<size_t>(old_encode_count) * stride);
    std::vector<void*> old_encode_work(old_encode_count);
    for (unsigned i = 0; i < old_encode_count; ++i)
        old_encode_work[i] = old_encode_storage.data() +
            static_cast<size_t>(i) * stride;
    require(leo_encode(stride, test.k, test.r, old_encode_count,
        &original_input[0], &old_encode_work[0]) == Leopard_Success,
        "legacy encode");

    leo2_codec_options codec_options = {};
    codec_options.struct_size = sizeof(codec_options);
    codec_options.flags = LEO2_CODEC_FORCE_SPECIALIZED_DECODE |
        ((case_index & 1u) != 0 ? LEO2_CODEC_FORCE_TILED_DECODE
                                : LEO2_CODEC_FORCE_MATERIALIZED_DECODE);
    leo2_codec* codec = NULL;
    require_result(leo2_codec_create(context, test.k, test.r,
        LEO2_PROFILE_LEGACY_HIGH_V1, test.field, &codec_options, &codec),
        "codec create");

    AlignedBytes new_recovery_storage(static_cast<size_t>(test.r) * stride);
    std::vector<void*> new_recovery_output(test.r);
    for (uint32_t i = 0; i < test.r; ++i)
        new_recovery_output[i] = new_recovery_storage.data() +
            static_cast<size_t>(i) * stride;
    size_t encode_scratch_bytes = 0;
    require_result(leo2_encode_scratch_size(
        codec, test.bytes, &encode_scratch_bytes), "encode scratch query");
    AlignedBytes encode_scratch(encode_scratch_bytes);
    require_result(leo2_encode(codec, test.bytes, &original_input[0],
        &new_recovery_output[0], encode_scratch.data(), encode_scratch.size()),
        "new encode");
    for (uint32_t i = 0; i < test.r; ++i)
    {
        require(memcmp(new_recovery_output[i], old_encode_work[i],
            test.bytes) == 0, "new parity differs from Leopard1");
        parity_bytes += test.bytes;
    }

    const uint32_t loss_count = std::min<uint32_t>(5,
        std::min(test.k, test.r));
    std::vector<uint32_t> missing(loss_count);
    for (uint32_t i = 0; i < loss_count; ++i)
        missing[i] = loss_count == 1 ? 0 :
            static_cast<uint32_t>((static_cast<uint64_t>(i) * (test.k - 1)) /
                (loss_count - 1));

    std::vector<const void*> old_original = original_input;
    std::vector<uint8_t> original_present(test.k, 1);
    for (uint32_t i = 0; i < loss_count; ++i)
    {
        old_original[missing[i]] = NULL;
        original_present[missing[i]] = 0;
    }
    std::vector<const void*> old_recovery(test.r);
    for (uint32_t i = 0; i < test.r; ++i)
        old_recovery[i] = old_encode_work[i];
    const unsigned old_decode_count = leo_decode_work_count(test.k, test.r);
    require(old_decode_count != 0, "legacy decode work count");
    AlignedBytes old_decode_storage(
        static_cast<size_t>(old_decode_count) * stride);
    std::vector<void*> old_decode_work(old_decode_count);
    for (unsigned i = 0; i < old_decode_count; ++i)
        old_decode_work[i] = old_decode_storage.data() +
            static_cast<size_t>(i) * stride;
    require(leo_decode(stride, test.k, test.r, old_decode_count,
        &old_original[0], &old_recovery[0], &old_decode_work[0]) ==
        Leopard_Success, "legacy decode");

    std::vector<uint8_t> recovery_present(test.r, 1);
    leo2_decode_plan* plan = NULL;
    require_result(leo2_decode_plan_create(codec, &original_present[0],
        &recovery_present[0], &plan), "decode plan create");
    size_t decode_scratch_bytes = 0;
    require_result(leo2_decode_plan_scratch_size(
        plan, test.bytes, &decode_scratch_bytes), "decode scratch query");
    AlignedBytes decode_scratch(decode_scratch_bytes);
    AlignedBytes restored_storage(static_cast<size_t>(test.k) * stride);
    std::vector<void*> restored(test.k, NULL);
    for (uint32_t i = 0; i < loss_count; ++i)
        restored[missing[i]] = restored_storage.data() +
            static_cast<size_t>(missing[i]) * stride;
    std::vector<const void*> new_recovery(test.r);
    for (uint32_t i = 0; i < test.r; ++i)
        new_recovery[i] = new_recovery_output[i];
#if defined(LEO2_ENABLE_TEST_HOOKS)
    const bool verify_materialized_cancellation =
        test.field == LEO2_FIELD_GF16 && (case_index & 1u) == 0;
    if (verify_materialized_cancellation)
        leopard::ff16::TestOnlyResetHighDecodeCounts();
#endif
    require_result(leo2_decode_plan_execute(plan, test.bytes,
        &old_original[0], &new_recovery[0], &restored[0],
        decode_scratch.data(), decode_scratch.size()), "new decode");
#if defined(LEO2_ENABLE_TEST_HOOKS)
    if (verify_materialized_cancellation)
    {
        const leopard::ff16::TestOnlyHighDecodeCounts stage_counts =
            leopard::ff16::TestOnlyGetHighDecodeCounts();
        require(stage_counts.receive_ifft_butterfly4_out_of_place == 0,
            "materialized GF16 unexpectedly fused receive rows");
        require(stage_counts.receive_copy_shards == test.k - loss_count,
            "materialized GF16 did not stage each later selected row once");
        require(stage_counts.syndrome_block_zero_ifft_elisions == 1,
            "materialized GF16 did not elide the block-zero inverse transform");
        require(stage_counts.syndrome_forward_transforms == 1 &&
                stage_counts.syndrome_forward_transform_elisions == 0,
            "materialized GF16 later syndrome forward count drifted");
        require(stage_counts.syndrome_block_zero_xor_shards == loss_count,
            "materialized GF16 block-zero raw-source XOR count drifted");
    }
#endif
    for (uint32_t i = 0; i < loss_count; ++i)
    {
        const uint32_t shard = missing[i];
        require(memcmp(restored[shard], old_decode_work[shard],
            test.bytes) == 0, "new decode differs from Leopard1");
        require(memcmp(restored[shard], original_input[shard],
            test.bytes) == 0, "restored original mismatch");
        restored_bytes += test.bytes;
    }

#if defined(LEO2_ENABLE_TEST_HOOKS)
    verify_gf8_one_loss_receive_stage(
        context, codec, test, original_input, new_recovery_output);
#endif

    leo2_decode_plan_destroy(plan);
    leo2_codec_destroy(codec);
}

} // namespace

int main()
{
    try
    {
        require(leo_init() == Leopard_Success, "Leopard initialization");
        leo2_context_options options = {};
        options.struct_size = sizeof(options);
        options.thread_count = 1;
        leo2_context* context = NULL;
        require_result(leo2_context_create(&options, &context),
            "context create");
        const Case cases[] = {
            { 3, 2, LEO2_FIELD_GF8, 65 },
            { 7, 3, LEO2_FIELD_GF8, 64 },
            { 17, 5, LEO2_FIELD_GF8, 129 },
            { 33, 9, LEO2_FIELD_GF8, 64 },
            { 65, 17, LEO2_FIELD_GF8, 257 },
            { 127, 31, LEO2_FIELD_GF8, 64 },
            { 191, 32, LEO2_FIELD_GF8, 129 },
            { 223, 16, LEO2_FIELD_GF8, 65 },
            { 240, 16, LEO2_FIELD_GF8, 64 },
            { 193, 65, LEO2_FIELD_GF16, 128 },
            { 257, 32, LEO2_FIELD_GF16, 128 },
            { 511, 64, LEO2_FIELD_GF16, 128 },
            { 1000, 200, LEO2_FIELD_GF16, 128 }
        };
        uint64_t parity_bytes = 0;
        uint64_t restored_bytes = 0;
        for (size_t i = 0; i < sizeof(cases) / sizeof(cases[0]); ++i)
            run_case(context, cases[i], i, parity_bytes, restored_bytes);
#if defined(LEO2_ENABLE_TEST_HOOKS)
        verify_block_zero_only_cancellation(context);
        verify_internal_cancellation_patterns();
#endif
        leo2_context_destroy(context);
        std::cout << "PASS high_pruned_legacy cases="
                  << sizeof(cases) / sizeof(cases[0])
                  << " parity_bytes=" << parity_bytes
                  << " restored_bytes=" << restored_bytes
                  << " workspaces=materialized,tiled" << std::endl;
        return 0;
    }
    catch (const std::exception& error)
    {
        std::cerr << "FAIL high_pruned_legacy: " << error.what() << std::endl;
        return 1;
    }
}
