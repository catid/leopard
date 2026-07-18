/*
    Copyright (c) 2017 Christopher A. Taylor.  All rights reserved.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the conditions in the Leopard-RS
    BSD license are met.  See the repository LICENSE file.
*/

#ifndef LEO2_ENABLE_TEST_HOOKS
#error "fuzz_pruned_transform requires LEO2_ENABLE_TEST_HOOKS"
#endif

#include "Leopard2Backend.h"
#include "Leopard2Plan.h"
#include "LeopardFF8.h"
#include "LeopardFF16.h"
#include "leopard.h"

#include <algorithm>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vector>

namespace {

void CheckImpl(bool condition, int line)
{
    if (!condition)
    {
        fprintf(stderr, "pruned transform fuzz invariant failed at line %d\n",
            line);
        abort();
    }
}

#define Check(condition) CheckImpl((condition), __LINE__)

uint64_t HashInput(const uint8_t* data, size_t size)
{
    uint64_t hash = UINT64_C(1469598103934665603);
    for (size_t i = 0; i < size; ++i)
    {
        hash ^= data[i];
        hash *= UINT64_C(1099511628211);
    }
    return hash;
}

class Random
{
public:
    explicit Random(uint64_t seed)
        : state_(seed ? seed : UINT64_C(0x9e3779b97f4a7c15))
    {}

    uint64_t Next()
    {
        uint64_t value = state_;
        value ^= value >> 12;
        value ^= value << 25;
        value ^= value >> 27;
        state_ = value;
        return value * UINT64_C(2685821657736338717);
    }

private:
    uint64_t state_;
};

struct Workspace
{
    std::vector<std::vector<uint8_t> > shards;
    std::vector<void*> pointers;

    Workspace(unsigned size, size_t bytes)
        : shards(size, std::vector<uint8_t>(bytes, 0))
        , pointers(size, NULL)
    {
        for (unsigned i = 0; i < size; ++i)
            pointers[i] = shards[i].data();
    }
};

unsigned PrefixExtent(const std::vector<uint8_t>& mask)
{
    unsigned extent = static_cast<unsigned>(mask.size());
    while (extent != 0 && mask[extent - 1] == 0)
        --extent;
    return extent;
}

std::vector<uint8_t> PackMask(const std::vector<uint8_t>& mask)
{
    std::vector<uint8_t> packed((mask.size() + 7u) / 8u, 0);
    for (size_t i = 0; i < mask.size(); ++i)
        if (mask[i])
            packed[i >> 3] |= static_cast<uint8_t>(1u << (i & 7u));
    return packed;
}

bool SamePlan(
    const leopard2_internal::PrunedTransformPlan& left,
    const leopard2_internal::PrunedTransformPlan& right)
{
    if (left.size != right.size || left.shift != right.shift ||
        left.zero_multiplier_log != right.zero_multiplier_log ||
        left.first_layer_multiplier_log !=
            right.first_layer_multiplier_log ||
        left.inverse != right.inverse ||
        left.input_mask != right.input_mask ||
        left.output_mask != right.output_mask ||
        left.fused_four_starts != right.fused_four_starts ||
        left.zero_outputs != right.zero_outputs ||
        left.inverse_source_prefix != right.inverse_source_prefix ||
        left.inverse_source_groups.size() !=
            right.inverse_source_groups.size() ||
        left.operations.size() != right.operations.size() ||
        left.full_butterfly_count != right.full_butterfly_count ||
        left.one_output_butterflies != right.one_output_butterflies ||
        left.input_zero_specializations != right.input_zero_specializations ||
        left.zero_multiplier_butterflies != right.zero_multiplier_butterflies ||
        left.one_multiplier_butterflies != right.one_multiplier_butterflies)
        return false;
    for (size_t i = 0; i < left.operations.size(); ++i)
    {
        const leopard2_internal::PrunedTransformOperation& a =
            left.operations[i];
        const leopard2_internal::PrunedTransformOperation& b =
            right.operations[i];
        if (a.x != b.x || a.y != b.y ||
            a.multiplier_log != b.multiplier_log || a.flags != b.flags)
            return false;
    }
    for (size_t i = 0; i < left.inverse_source_groups.size(); ++i)
    {
        const leopard2_internal::PrunedInverseSourceGroup& a =
            left.inverse_source_groups[i];
        const leopard2_internal::PrunedInverseSourceGroup& b =
            right.inverse_source_groups[i];
        if (a.multiplier_log01 != b.multiplier_log01 ||
            a.multiplier_log23 != b.multiplier_log23 ||
            a.multiplier_log02 != b.multiplier_log02)
            return false;
    }
    return true;
}

#ifdef LEO_HAS_FF8
struct GF8
{
    static unsigned Order() { return 256; }
    static bool Prepare(
        unsigned size, unsigned shift, bool inverse,
        const uint8_t* input, const uint8_t* output,
        leopard2_internal::PrunedTransformPlan& plan)
    {
        return leopard::ff8::PreparePrunedTransformPlan(
            size, shift, inverse, input, output, plan);
    }
    static bool PrepareSparse(
        unsigned size, unsigned shift, uint8_t* dependency,
        size_t dependency_bytes, uint8_t* operations, size_t operation_bytes,
        leopard2_internal::SparseForwardPlanStats& stats)
    {
        return leopard::ff8::PrepareSparseForwardPlan(
            size, shift, dependency, dependency_bytes,
            operations, operation_bytes, stats);
    }
    static bool ExecuteSparse(
        const leopard::backend::Ops& ops, uint64_t bytes,
        unsigned size, unsigned shift, const uint8_t* operations,
        size_t operation_bytes, void** work)
    {
        return leopard::ff8::ExecuteSparseForwardPlan(
            ops, bytes, size, shift, operations, operation_bytes, work);
    }
    static bool ExecuteSparseSources(
        const leopard::backend::Ops& ops, uint64_t bytes,
        unsigned size, unsigned shift, const uint8_t* operations,
        size_t operation_bytes, void* const* source, void** work)
    {
        return leopard::ff8::ExecuteSparseForwardPlanFromSources(
            ops, bytes, size, shift, operations, operation_bytes,
            source, work);
    }
    static void Full(
        const leopard::backend::Ops& ops,
        bool inverse,
        uint64_t bytes,
        unsigned size,
        unsigned shift,
        unsigned input_prefix,
        void** work)
    {
        if (inverse)
            leopard::ff8::TestOnlyLchInverse(
                ops, bytes, size, shift, input_prefix, work);
        else
            leopard::ff8::TestOnlyLchForward(
                ops, bytes, size, shift, size, work);
    }
};
#endif

#ifdef LEO_HAS_FF16
struct GF16
{
    static unsigned Order() { return 65536; }
    static bool Prepare(
        unsigned size, unsigned shift, bool inverse,
        const uint8_t* input, const uint8_t* output,
        leopard2_internal::PrunedTransformPlan& plan)
    {
        return leopard::ff16::PreparePrunedTransformPlan(
            size, shift, inverse, input, output, plan);
    }
    static bool PrepareSparse(
        unsigned size, unsigned shift, uint8_t* dependency,
        size_t dependency_bytes, uint8_t* operations, size_t operation_bytes,
        leopard2_internal::SparseForwardPlanStats& stats)
    {
        return leopard::ff16::PrepareSparseForwardPlan(
            size, shift, dependency, dependency_bytes,
            operations, operation_bytes, stats);
    }
    static bool ExecuteSparse(
        const leopard::backend::Ops& ops, uint64_t bytes,
        unsigned size, unsigned shift, const uint8_t* operations,
        size_t operation_bytes, void** work)
    {
        return leopard::ff16::ExecuteSparseForwardPlan(
            ops, bytes, size, shift, operations, operation_bytes, work);
    }
    static bool ExecuteSparseSources(
        const leopard::backend::Ops& ops, uint64_t bytes,
        unsigned size, unsigned shift, const uint8_t* operations,
        size_t operation_bytes, void* const* source, void** work)
    {
        return leopard::ff16::ExecuteSparseForwardPlanFromSources(
            ops, bytes, size, shift, operations, operation_bytes,
            source, work);
    }
    static void Full(
        const leopard::backend::Ops& ops,
        bool inverse,
        uint64_t bytes,
        unsigned size,
        unsigned shift,
        unsigned input_prefix,
        void** work)
    {
        if (inverse)
            leopard::ff16::TestOnlyLchInverse(
                ops, bytes, size, shift, input_prefix, work);
        else
            leopard::ff16::TestOnlyLchForward(
                ops, bytes, size, shift, size, work);
    }
};
#endif

template<class Field>
void RunCase(
    const leopard::backend::Ops& ops,
    Random& random,
    const uint8_t* data,
    size_t data_size,
    bool inverse,
    unsigned size,
    unsigned shift,
    size_t bytes)
{
    std::vector<uint8_t> input_mask(size, 0);
    std::vector<uint8_t> output_mask(size, 0);
    const unsigned density = 1u + (data_size > 4 ? data[4] % 7u : 3u);
    for (unsigned i = 0; i < size; ++i)
    {
        input_mask[i] = static_cast<uint8_t>(random.Next() % density == 0);
        output_mask[i] = static_cast<uint8_t>(random.Next() % density == 0);
    }
    output_mask[random.Next() % size] = 1;

    leopard2_internal::PrunedTransformPlan fused_plan;
    Check(Field::Prepare(
        size, shift, inverse, input_mask.data(), output_mask.data(),
        fused_plan));
    const leopard2_internal::PrunedTransformPlan plan_snapshot(fused_plan);

    // A failed setup must not partially publish over a reusable immutable
    // plan.  Exercise every malformed byte class through libFuzzer's input.
    const uint8_t malformed_value = static_cast<uint8_t>(
        2u + (data_size > 5 ? data[5] % 254u : 0u));
    std::vector<uint8_t> malformed(input_mask);
    malformed[random.Next() % size] = malformed_value;
    Check(!Field::Prepare(
        size, shift, inverse, malformed.data(), output_mask.data(),
        fused_plan));
    Check(SamePlan(fused_plan, plan_snapshot));
    malformed = output_mask;
    malformed[random.Next() % size] = malformed_value;
    Check(!Field::Prepare(
        size, shift, inverse, input_mask.data(), malformed.data(),
        fused_plan));
    Check(SamePlan(fused_plan, plan_snapshot));

    leopard2_internal::PrunedTransformPlan flat_plan(fused_plan);
    std::vector<uint32_t>().swap(flat_plan.fused_four_starts);

    Workspace initial(size, bytes);
    for (unsigned shard = 0; shard < size; ++shard)
    {
        if (!input_mask[shard])
            continue;
        for (size_t offset = 0; offset < bytes; ++offset)
            initial.shards[shard][offset] =
                static_cast<uint8_t>(random.Next());
    }
    Workspace full(initial);
    Workspace flat(initial);
    Workspace fused(initial);
    // The implicitly copied pointer arrays still refer to the source object.
    // Rebind them before entering any transform.
    for (unsigned i = 0; i < size; ++i)
    {
        full.pointers[i] = full.shards[i].data();
        flat.pointers[i] = flat.shards[i].data();
        fused.pointers[i] = fused.shards[i].data();
    }

    Field::Full(
        ops, inverse, bytes, size, shift, PrefixExtent(input_mask),
        full.pointers.data());
    Check(leopard2_internal::ExecutePrunedTransformPlan(
        ops, bytes, flat_plan, flat.pointers.data()));
    Check(leopard2_internal::ExecutePrunedTransformPlan(
        ops, bytes, fused_plan, fused.pointers.data()));
    for (unsigned i = 0; i < size; ++i)
    {
        if (output_mask[i])
        {
            Check(full.shards[i] == flat.shards[i]);
            Check(full.shards[i] == fused.shards[i]);
        }
    }

    if (!inverse)
    {
        std::vector<uint8_t> complete_input(size, 1);
        leopard2_internal::PrunedTransformPlan source_plan;
        Check(Field::Prepare(
            size, shift, false, complete_input.data(), output_mask.data(),
            source_plan));
        Workspace source(size, bytes);
        for (unsigned shard = 0; shard < size; ++shard)
            for (size_t offset = 0; offset < bytes; ++offset)
                source.shards[shard][offset] =
                    static_cast<uint8_t>(random.Next());
        Workspace source_snapshot(source);
        Workspace source_expected(source);
        Workspace source_actual(size, bytes);
        for (unsigned i = 0; i < size; ++i)
        {
            source_snapshot.pointers[i] = source_snapshot.shards[i].data();
            source_expected.pointers[i] = source_expected.shards[i].data();
        }
        Field::Full(
            ops, false, bytes, size, shift, size,
            source_expected.pointers.data());
        Check(leopard2_internal::ExecutePrunedForwardTransformPlanFromSources(
            ops, bytes, source_plan, source.pointers.data(),
            source_actual.pointers.data()));
        Check(source.shards == source_snapshot.shards);
        for (unsigned i = 0; i < size; ++i)
            if (output_mask[i])
                Check(source_actual.shards[i] == source_expected.shards[i]);

        std::vector<uint8_t> dependency = PackMask(output_mask);
        std::vector<uint8_t> operation_masks(
            leopard2_internal::SparseForwardRetainedBytes(size), 0xa5);
        leopard2_internal::SparseForwardPlanStats sparse_stats;
        Check(Field::PrepareSparse(
            size, shift, dependency.data(), dependency.size(),
            operation_masks.data(), operation_masks.size(), sparse_stats));
        Check(sparse_stats.retained_butterfly_count ==
            leopard2_internal::CountSparseForwardRetainedButterflies(
                size, operation_masks.data(), operation_masks.size()));

        Workspace sparse_in_place(source);
        Workspace sparse_source(size, bytes);
        for (unsigned i = 0; i < size; ++i)
            sparse_in_place.pointers[i] = sparse_in_place.shards[i].data();
        Check(Field::ExecuteSparse(
            ops, bytes, size, shift, operation_masks.data(),
            operation_masks.size(), sparse_in_place.pointers.data()));
        Check(Field::ExecuteSparseSources(
            ops, bytes, size, shift, operation_masks.data(),
            operation_masks.size(), source.pointers.data(),
            sparse_source.pointers.data()));
        Check(source.shards == source_snapshot.shards);
        for (unsigned i = 0; i < size; ++i)
        {
            if (output_mask[i])
            {
                Check(sparse_in_place.shards[i] == source_expected.shards[i]);
                Check(sparse_source.shards[i] == source_expected.shards[i]);
            }
        }

        leopard2_internal::PrunedTransformPlan broken_source(source_plan);
        broken_source.inverse = true;
        Check(!leopard2_internal::ExecutePrunedForwardTransformPlanFromSources(
            ops, bytes, broken_source, source.pointers.data(),
            source_actual.pointers.data()));
        broken_source = source_plan;
        broken_source.input_mask[random.Next() % size] = 0;
        Check(!leopard2_internal::ExecutePrunedForwardTransformPlanFromSources(
            ops, bytes, broken_source, source.pointers.data(),
            source_actual.pointers.data()));
    }
    else
    {
        const unsigned prefix = 1U +
            static_cast<unsigned>(random.Next() % (size - 1U));
        std::vector<uint8_t> prefix_input(size, 0);
        std::vector<uint8_t> complete_output(size, 1);
        std::fill(prefix_input.begin(), prefix_input.begin() + prefix, 1);
        leopard2_internal::PrunedTransformPlan source_plan;
        Check(Field::Prepare(
            size, shift, true, prefix_input.data(), complete_output.data(),
            source_plan));
        Workspace source(size, bytes);
        for (unsigned shard = 0; shard < prefix; ++shard)
            for (size_t offset = 0; offset < bytes; ++offset)
                source.shards[shard][offset] =
                    static_cast<uint8_t>(random.Next());
        Workspace source_snapshot(source);
        Workspace source_expected(source);
        Workspace source_actual(size, bytes);
        for (unsigned i = 0; i < size; ++i)
        {
            source_snapshot.pointers[i] = source_snapshot.shards[i].data();
            source_expected.pointers[i] = source_expected.shards[i].data();
        }
        Field::Full(
            ops, true, bytes, size, shift, prefix,
            source_expected.pointers.data());
        std::vector<void*> prefix_pointers(prefix, NULL);
        for (unsigned i = 0; i < prefix; ++i)
            prefix_pointers[i] = source.pointers[i];
        Check(leopard2_internal::ExecutePrunedInverseTransformPlanFromSources(
            ops, bytes, source_plan, prefix_pointers.data(),
            source_actual.pointers.data()));
        Check(source.shards == source_snapshot.shards);
        Check(source_actual.shards == source_expected.shards);

        leopard2_internal::PrunedTransformPlan broken_source(source_plan);
        broken_source.output_mask[random.Next() % size] = 0;
        Check(!leopard2_internal::ExecutePrunedInverseTransformPlanFromSources(
            ops, bytes, broken_source, prefix_pointers.data(),
            source_actual.pointers.data()));
        broken_source = source_plan;
        broken_source.inverse_source_prefix = size;
        Check(!leopard2_internal::ExecutePrunedInverseTransformPlanFromSources(
            ops, bytes, broken_source, prefix_pointers.data(),
            source_actual.pointers.data()));
    }

    // The executor validates every stored index and multiplier before using
    // it.  Mutations must fail rather than touch a shard outside the plan.
    leopard2_internal::PrunedTransformPlan broken(fused_plan);
    if (!broken.operations.empty())
    {
        broken.operations[0].x = size;
        Check(!leopard2_internal::ExecutePrunedTransformPlan(
            ops, bytes, broken, fused.pointers.data()));
        if (broken.zero_multiplier_log != UINT16_MAX)
        {
            broken = fused_plan;
            broken.operations[0].multiplier_log = static_cast<uint16_t>(
                broken.zero_multiplier_log + 1u);
            Check(!leopard2_internal::ExecutePrunedTransformPlan(
                ops, bytes, broken, fused.pointers.data()));
        }
    }
    broken = fused_plan;
    broken.zero_outputs.push_back(size);
    Check(!leopard2_internal::ExecutePrunedTransformPlan(
        ops, bytes, broken, fused.pointers.data()));
    broken = fused_plan;
    broken.fused_four_starts.insert(
        broken.fused_four_starts.begin(),
        static_cast<uint32_t>(broken.operations.size()));
    Check(!leopard2_internal::ExecutePrunedTransformPlan(
        ops, bytes, broken, fused.pointers.data()));
}

} // namespace

extern "C" int LLVMFuzzerTestOneInput(const uint8_t* data, size_t size)
{
    if (!data || size < 6)
        return 0;
    try
    {
        Check(leo_init() == Leopard_Success);
        Random random(HashInput(data, size));
#if defined(LEO_HAS_FF8) && defined(LEO_HAS_FF16)
        const bool gf16 = 0 != (data[0] & 1u);
#elif defined(LEO_HAS_FF16)
        const bool gf16 = true;
#else
        const bool gf16 = false;
#endif
        const unsigned log_size = 1u + data[1] % 8u;
        const unsigned transform_size = 1u << log_size;
#if defined(LEO_HAS_FF8) && defined(LEO_HAS_FF16)
        const unsigned order = gf16 ? GF16::Order() : GF8::Order();
#elif defined(LEO_HAS_FF16)
        const unsigned order = GF16::Order();
#else
        const unsigned order = GF8::Order();
#endif
        const unsigned block_count = order / transform_size;
        const unsigned shift = static_cast<unsigned>(
            random.Next() % block_count) * transform_size;
        const bool inverse = 0 != (data[2] & 1u);
        size_t bytes = 1u +
            (static_cast<size_t>(data[3]) |
             (static_cast<size_t>(data[4]) << 8)) % 257u;
        if (gf16)
            bytes = (bytes + 1u) & ~static_cast<size_t>(1u);

        const leo2_backend requested[] = {
            LEO2_BACKEND_SCALAR,
            LEO2_BACKEND_SSSE3,
            LEO2_BACKEND_AVX2
        };
        leopard::backend::QualificationStatus status =
            leopard::backend::QualificationAvailable;
        const leopard::backend::Ops* ops = leopard::backend::GetQualifiedOps(
            requested[data[5] % 3u], &status);
        if (!ops)
        {
            ops = leopard::backend::GetQualifiedOps(
                LEO2_BACKEND_SCALAR, &status);
            Check(ops != NULL);
        }

#if defined(LEO_HAS_FF8) && defined(LEO_HAS_FF16)
        if (gf16)
            RunCase<GF16>(
                *ops, random, data, size, inverse, transform_size,
                shift, bytes);
        else
            RunCase<GF8>(
                *ops, random, data, size, inverse, transform_size,
                shift, bytes);
#elif defined(LEO_HAS_FF16)
        RunCase<GF16>(
            *ops, random, data, size, inverse, transform_size,
            shift, bytes);
#else
        RunCase<GF8>(
            *ops, random, data, size, inverse, transform_size,
            shift, bytes);
#endif
    }
    catch (...)
    {
        abort();
    }
    return 0;
}
