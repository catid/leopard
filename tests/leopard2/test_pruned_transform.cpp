/*
    Copyright (c) 2017 Christopher A. Taylor.  All rights reserved.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the conditions in the Leopard-RS
    BSD license are met.  See the repository LICENSE file.
*/

#include "Leopard2Backend.h"
#include "Leopard2Plan.h"
#include "LeopardFF8.h"
#include "LeopardFF16.h"
#include "leopard.h"
#include "direct_oracle.h"

#include <atomic>
#include <cstring>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <thread>
#include <vector>

namespace {

struct TestCounts
{
    uint64_t plans;
    uint64_t executions;
    uint64_t compared_bytes;
    uint64_t direct_symbols;
    uint64_t backends;
    uint64_t fused_four_descriptors;
    uint64_t execution_steps;
    uint64_t max_plan_bytes;

    TestCounts()
        : plans(0), executions(0), compared_bytes(0), direct_symbols(0)
        , backends(0), fused_four_descriptors(0), execution_steps(0)
        , max_plan_bytes(0)
    {}
};

void require(bool condition, const std::string& message)
{
    if (!condition)
        throw std::runtime_error(message);
}

uint64_t plan_storage_bytes(
    const leopard2_internal::PrunedTransformPlan& plan)
{
    return sizeof(plan) +
        plan.input_mask.capacity() * sizeof(plan.input_mask[0]) +
        plan.output_mask.capacity() * sizeof(plan.output_mask[0]) +
        plan.operations.capacity() * sizeof(plan.operations[0]) +
        plan.fused_four_starts.capacity() *
            sizeof(plan.fused_four_starts[0]) +
        plan.zero_outputs.capacity() * sizeof(plan.zero_outputs[0]);
}

uint64_t mix64(uint64_t value)
{
    value ^= value >> 30;
    value *= UINT64_C(0xbf58476d1ce4e5b9);
    value ^= value >> 27;
    value *= UINT64_C(0x94d049bb133111eb);
    return value ^ (value >> 31);
}

std::vector<void*> pointers(std::vector<std::vector<uint8_t> >& shards)
{
    std::vector<void*> result(shards.size(), NULL);
    for (size_t i = 0; i < shards.size(); ++i)
        result[i] = shards[i].data();
    return result;
}

std::vector<void*> source_pointers(
    std::vector<std::vector<uint8_t> >& shards)
{
    // The internal executor accepts mutable pointer values for compatibility
    // with the existing transform plumbing, but treats every pointed-to byte
    // as immutable.
    return pointers(shards);
}

std::vector<std::vector<uint8_t> > make_input(
    unsigned size,
    uint64_t bytes,
    const std::vector<uint8_t>& active,
    uint64_t seed)
{
    std::vector<std::vector<uint8_t> > shards(
        size, std::vector<uint8_t>(static_cast<size_t>(bytes), 0));
    for (unsigned shard = 0; shard < size; ++shard)
    {
        if (!active[shard])
            continue;
        for (uint64_t offset = 0; offset < bytes; ++offset)
        {
            const uint64_t value = mix64(
                seed ^ (static_cast<uint64_t>(shard) << 32) ^ offset);
            shards[shard][static_cast<size_t>(offset)] =
                static_cast<uint8_t>(value);
        }
    }
    return shards;
}

std::vector<std::vector<uint8_t> > make_masks(unsigned size)
{
    std::vector<std::vector<uint8_t> > masks;
    masks.push_back(std::vector<uint8_t>(size, 0));

    std::vector<uint8_t> first(size, 0);
    first[0] = 1;
    masks.push_back(first);

    std::vector<uint8_t> prefix(size, 0);
    const unsigned prefix_count = size / 2 + 1;
    for (unsigned i = 0; i < prefix_count; ++i)
        prefix[i] = 1;
    masks.push_back(prefix);

    std::vector<uint8_t> sparse(size, 0);
    for (unsigned i = 0; i < size; ++i)
        sparse[i] = static_cast<uint8_t>(((i * 5U + 3U) % 11U) < 3U);
    sparse[size - 1] = 1;
    masks.push_back(sparse);

    masks.push_back(std::vector<uint8_t>(size, 1));
    return masks;
}

std::vector<unsigned> shifts(unsigned order, unsigned size)
{
    std::vector<unsigned> result;
    result.push_back(0);
    if (size < order)
        result.push_back(size);
    if (order - size != 0 && order - size != size)
        result.push_back(order - size);
    return result;
}

unsigned ceil_power_of_two(unsigned value)
{
    unsigned result = 1;
    while (result < value)
        result <<= 1;
    return result;
}

unsigned exact_log2(unsigned value)
{
    unsigned result = 0;
    while ((static_cast<unsigned>(1) << result) < value)
        ++result;
    require((static_cast<unsigned>(1) << result) == value,
        "non-dyadic direct-oracle size");
    return result;
}

struct GF8
{
    static unsigned order() { return leopard::ff8::kOrder; }
    static const leopard2_test::BinaryField& direct_field()
    {
        static const leopard2_test::BinaryField field =
            leopard2_test::make_legacy_gf8();
        return field;
    }
    static uint64_t symbol_bytes() { return 1; }
    static void set_symbol(
        std::vector<uint8_t>& shard, leopard2_test::Element value)
    {
        shard[0] = static_cast<uint8_t>(value);
    }
    static leopard2_test::Element get_symbol(
        const std::vector<uint8_t>& shard)
    {
        return shard[0];
    }
    static bool prepare(
        unsigned size, unsigned shift, bool inverse,
        const uint8_t* input, const uint8_t* output,
        leopard2_internal::PrunedTransformPlan& plan)
    {
        return leopard::ff8::PreparePrunedTransformPlan(
            size, shift, inverse, input, output, plan);
    }
    static void full(
        const leopard::backend::Ops& ops, bool inverse, uint64_t bytes,
        unsigned size, unsigned shift, void** work)
    {
        if (inverse)
            leopard::ff8::TestOnlyLchInverse(
                ops, bytes, size, shift, size, work);
        else
            leopard::ff8::TestOnlyLchForward(
                ops, bytes, size, shift, size, work);
    }
};

struct GF16
{
    static unsigned order() { return leopard::ff16::kOrder; }
    static const leopard2_test::BinaryField& direct_field()
    {
        static const leopard2_test::BinaryField field =
            leopard2_test::make_legacy_gf16();
        return field;
    }
    static uint64_t symbol_bytes() { return 2; }
    static void set_symbol(
        std::vector<uint8_t>& shard, leopard2_test::Element value)
    {
        shard[0] = static_cast<uint8_t>(value);
        shard[1] = static_cast<uint8_t>(value >> 8);
    }
    static leopard2_test::Element get_symbol(
        const std::vector<uint8_t>& shard)
    {
        return static_cast<leopard2_test::Element>(shard[0] |
            (static_cast<unsigned>(shard[1]) << 8));
    }
    static bool prepare(
        unsigned size, unsigned shift, bool inverse,
        const uint8_t* input, const uint8_t* output,
        leopard2_internal::PrunedTransformPlan& plan)
    {
        return leopard::ff16::PreparePrunedTransformPlan(
            size, shift, inverse, input, output, plan);
    }
    static void full(
        const leopard::backend::Ops& ops, bool inverse, uint64_t bytes,
        unsigned size, unsigned shift, void** work)
    {
        if (inverse)
            leopard::ff16::TestOnlyLchInverse(
                ops, bytes, size, shift, size, work);
        else
            leopard::ff16::TestOnlyLchForward(
                ops, bytes, size, shift, size, work);
    }
};

template<class Field>
void run_case(
    const leopard::backend::Ops& ops,
    unsigned size,
    unsigned shift,
    bool inverse,
    uint64_t bytes,
    const std::vector<uint8_t>& input_mask,
    const std::vector<uint8_t>& output_mask,
    uint64_t seed,
    TestCounts& counts)
{
    leopard2_internal::PrunedTransformPlan plan;
    require(Field::prepare(
            size, shift, inverse, input_mask.data(), output_mask.data(), plan),
        "pruned plan construction failed");
    require(plan.size == size && plan.shift == shift &&
            plan.inverse == inverse,
        "pruned plan identity mismatch");
    require(plan.operations.size() <= plan.full_butterfly_count,
        "pruned plan exceeds padded transform");
    counts.fused_four_descriptors += plan.fused_four_starts.size();
    counts.execution_steps += plan.operations.size() -
        plan.fused_four_starts.size() * 3U;
    const uint64_t plan_bytes = plan_storage_bytes(plan);
    if (plan_bytes > counts.max_plan_bytes)
        counts.max_plan_bytes = plan_bytes;

    std::vector<std::vector<uint8_t> > initial =
        make_input(size, bytes, input_mask, seed);
    std::vector<std::vector<uint8_t> > expected(initial);
    std::vector<std::vector<uint8_t> > actual(initial);
    std::vector<void*> expected_pointers = pointers(expected);
    std::vector<void*> actual_pointers = pointers(actual);
    Field::full(
        ops, inverse, bytes, size, shift, expected_pointers.data());
    require(leopard2_internal::ExecutePrunedTransformPlan(
            ops, bytes, plan, actual_pointers.data()),
        "pruned transform execution failed");

    for (unsigned coordinate = 0; coordinate < size; ++coordinate)
    {
        if (!output_mask[coordinate])
            continue;
        if (std::memcmp(actual[coordinate].data(),
                expected[coordinate].data(), static_cast<size_t>(bytes)) != 0)
        {
            std::ostringstream stream;
            stream << "pruned transform mismatch backend=" << ops.name
                   << " size=" << size << " shift=" << shift
                   << " inverse=" << inverse << " bytes=" << bytes
                   << " coordinate=" << coordinate;
            throw std::runtime_error(stream.str());
        }
        counts.compared_bytes += bytes;
    }
    ++counts.plans;
    ++counts.executions;

    bool complete_input = !inverse;
    for (unsigned coordinate = 0; coordinate < size; ++coordinate)
        complete_input = complete_input && input_mask[coordinate] != 0;
    if (complete_input)
    {
        const std::vector<std::vector<uint8_t> > source_snapshot(initial);
        std::vector<std::vector<uint8_t> > from_source(
            size, std::vector<uint8_t>(static_cast<size_t>(bytes), 0xa5));
        std::vector<void*> immutable_source = source_pointers(initial);
        std::vector<void*> from_source_pointers = pointers(from_source);
        require(
            leopard2_internal::ExecutePrunedForwardTransformPlanFromSources(
                ops, bytes, plan, immutable_source.data(),
                from_source_pointers.data()),
            "immutable-source pruned transform execution failed");
        require(initial == source_snapshot,
            "immutable-source pruned transform changed its coefficients");
        for (unsigned coordinate = 0; coordinate < size; ++coordinate)
        {
            if (!output_mask[coordinate])
                continue;
            require(from_source[coordinate] == expected[coordinate],
                "immutable-source pruned transform output mismatch");
            counts.compared_bytes += bytes;
        }
        ++counts.executions;
    }
}

template<class Field>
void run_direct_oracle_case(
    const leopard::backend::Ops& ops,
    unsigned size,
    unsigned shift,
    bool inverse,
    const std::vector<uint8_t>& input_mask,
    const std::vector<uint8_t>& output_mask,
    uint64_t seed,
    TestCounts& counts)
{
    typedef leopard2_test::Element Element;
    const leopard2_test::BinaryField& field = Field::direct_field();
    std::vector<Element> input(size, 0);
    std::vector<std::vector<uint8_t> > actual(
        size, std::vector<uint8_t>(
            static_cast<size_t>(Field::symbol_bytes()), 0));
    for (unsigned i = 0; i < size; ++i)
    {
        if (!input_mask[i])
            continue;
        input[i] = static_cast<Element>(
            mix64(seed ^ (static_cast<uint64_t>(i) << 24)) &
            (Field::order() - 1U));
        Field::set_symbol(actual[i], input[i]);
    }

    std::vector<Element> expected(size, 0);
    const leopard2_test::LchBasis basis =
        leopard2_test::make_lch_basis(field, exact_log2(size));
    if (inverse)
    {
        std::vector<Element> points(size, 0);
        for (unsigned i = 0; i < size; ++i)
            points[i] = static_cast<Element>(shift ^ i);
        const leopard2_test::Polynomial polynomial =
            leopard2_test::lagrange_interpolate(field, points, input);
        expected = leopard2_test::polynomial_to_lch_coefficients(
            field, basis, polynomial);
    }
    else
    {
        const leopard2_test::Polynomial polynomial =
            leopard2_test::lch_coefficients_to_polynomial(
                field, basis, input);
        for (unsigned i = 0; i < size; ++i)
            expected[i] = leopard2_test::polynomial_evaluate(
                field, polynomial, static_cast<Element>(shift ^ i));
    }

    leopard2_internal::PrunedTransformPlan plan;
    require(Field::prepare(
            size, shift, inverse, input_mask.data(), output_mask.data(), plan),
        "direct-oracle plan construction failed");
    std::vector<void*> actual_pointers = pointers(actual);
    require(leopard2_internal::ExecutePrunedTransformPlan(
            ops, Field::symbol_bytes(), plan, actual_pointers.data()),
        "direct-oracle plan execution failed");
    for (unsigned i = 0; i < size; ++i)
    {
        if (!output_mask[i])
            continue;
        if (Field::get_symbol(actual[i]) != expected[i])
        {
            std::ostringstream stream;
            stream << "pruned transform differs from direct polynomial"
                   << " backend=" << ops.name << " size=" << size
                   << " shift=" << shift << " inverse=" << inverse
                   << " coordinate=" << i;
            throw std::runtime_error(stream.str());
        }
        ++counts.direct_symbols;
        counts.compared_bytes += Field::symbol_bytes();
    }
    ++counts.plans;
    ++counts.executions;
}

template<class Field>
void test_matrix(
    const leopard::backend::Ops& ops,
    const unsigned* sizes,
    size_t size_count,
    const uint64_t* byte_counts,
    size_t byte_count,
    TestCounts& counts)
{
    for (size_t s = 0; s < size_count; ++s)
    {
        const unsigned size = sizes[s];
        const std::vector<std::vector<uint8_t> > masks = make_masks(size);
        const std::vector<unsigned> transform_shifts =
            shifts(Field::order(), size);
        for (size_t q = 0; q < transform_shifts.size(); ++q)
        {
            for (unsigned inverse = 0; inverse < 2; ++inverse)
            {
                for (size_t b = 0; b < byte_count; ++b)
                {
                    // Diagonal and crossed mask pairs cover empty, prefix,
                    // sparse, and dense dependencies without an excessive
                    // Cartesian test matrix.
                    for (size_t pattern = 0; pattern < masks.size(); ++pattern)
                    {
                        const size_t output_pattern =
                            (pattern * 3U + inverse + q) % masks.size();
                        run_case<Field>(
                            ops, size, transform_shifts[q], inverse != 0,
                            byte_counts[b], masks[pattern],
                            masks[output_pattern],
                            UINT64_C(0x4331425052554e45) ^
                                (static_cast<uint64_t>(size) << 40) ^
                                (static_cast<uint64_t>(q) << 32) ^
                                (static_cast<uint64_t>(b) << 24) ^
                                (static_cast<uint64_t>(pattern) << 8) ^
                                inverse,
                            counts);
                    }
                }
            }
        }
    }
}

template<class Field>
void test_exhaustive_small_masks(
    const leopard::backend::Ops& ops,
    uint64_t bytes,
    TestCounts& counts)
{
    const unsigned sizes[] = { 2, 4 };
    for (size_t s = 0; s < sizeof(sizes) / sizeof(sizes[0]); ++s)
    {
        const unsigned size = sizes[s];
        const uint32_t mask_count = static_cast<uint32_t>(1U) << size;
        const unsigned test_shifts[] = { 0, Field::order() - size };
        for (size_t q = 0;
             q < sizeof(test_shifts) / sizeof(test_shifts[0]); ++q)
        {
            for (unsigned inverse = 0; inverse < 2; ++inverse)
            {
                for (uint32_t input_bits = 0;
                     input_bits < mask_count; ++input_bits)
                {
                    std::vector<uint8_t> input(size, 0);
                    for (unsigned i = 0; i < size; ++i)
                        input[i] = static_cast<uint8_t>(
                            (input_bits >> i) & 1U);
                    for (uint32_t output_bits = 0;
                         output_bits < mask_count; ++output_bits)
                    {
                        std::vector<uint8_t> output(size, 0);
                        for (unsigned i = 0; i < size; ++i)
                            output[i] = static_cast<uint8_t>(
                                (output_bits >> i) & 1U);
                        run_case<Field>(
                            ops, size, test_shifts[q], inverse != 0, bytes,
                            input, output,
                            UINT64_C(0x4558484155535449) ^
                                (static_cast<uint64_t>(size) << 48) ^
                                (static_cast<uint64_t>(q) << 40) ^
                                (static_cast<uint64_t>(input_bits) << 20) ^
                                (static_cast<uint64_t>(output_bits) << 4) ^
                                inverse,
                            counts);
                    }
                }
            }
        }
    }
}

template<class Field>
void test_direct_oracle_masks(
    const leopard::backend::Ops& ops,
    TestCounts& counts)
{
    const unsigned sizes[] = { 2, 4 };
    for (size_t s = 0; s < sizeof(sizes) / sizeof(sizes[0]); ++s)
    {
        const unsigned size = sizes[s];
        const uint32_t mask_count = static_cast<uint32_t>(1U) << size;
        const unsigned test_shifts[] = { 0, Field::order() - size };
        for (size_t q = 0;
             q < sizeof(test_shifts) / sizeof(test_shifts[0]); ++q)
            for (unsigned inverse = 0; inverse < 2; ++inverse)
                for (uint32_t input_bits = 0;
                     input_bits < mask_count; ++input_bits)
                    for (uint32_t output_bits = 0;
                         output_bits < mask_count; ++output_bits)
                    {
                        std::vector<uint8_t> input(size, 0);
                        std::vector<uint8_t> output(size, 0);
                        for (unsigned i = 0; i < size; ++i)
                        {
                            input[i] = static_cast<uint8_t>(
                                (input_bits >> i) & 1U);
                            output[i] = static_cast<uint8_t>(
                                (output_bits >> i) & 1U);
                        }
                        run_direct_oracle_case<Field>(
                            ops, size, test_shifts[q], inverse != 0,
                            input, output,
                            UINT64_C(0x4449524543544c43) ^
                                (static_cast<uint64_t>(size) << 48) ^
                                (static_cast<uint64_t>(q) << 40) ^
                                (static_cast<uint64_t>(input_bits) << 20) ^
                                (static_cast<uint64_t>(output_bits) << 4) ^
                                inverse,
                            counts);
                    }
    }

    const unsigned size = 8;
    const std::vector<std::vector<uint8_t> > masks = make_masks(size);
    const unsigned test_shifts[] = { 0, Field::order() - size };
    for (size_t q = 0;
         q < sizeof(test_shifts) / sizeof(test_shifts[0]); ++q)
        for (unsigned inverse = 0; inverse < 2; ++inverse)
            for (size_t pattern = 0; pattern < masks.size(); ++pattern)
                run_direct_oracle_case<Field>(
                    ops, size, test_shifts[q], inverse != 0,
                    masks[pattern], masks[(pattern * 2U + q + inverse) %
                        masks.size()],
                    UINT64_C(0x4449524543544e38) ^
                        (static_cast<uint64_t>(q) << 20) ^
                        (static_cast<uint64_t>(pattern) << 4) ^ inverse,
                    counts);
}

template<class Field>
void test_profile_masks(
    const leopard::backend::Ops& ops,
    unsigned k,
    unsigned r,
    uint64_t bytes,
    TestCounts& counts)
{
    const unsigned t = ceil_power_of_two(r);
    const unsigned high_parent = ceil_power_of_two(k + t);
    require(high_parent <= Field::order(), "high profile exceeds field");
    std::vector<uint8_t> all_t(t, 1);
    std::vector<uint8_t> high_partial(t, 0);
    unsigned final_message_count = k % t;
    if (final_message_count == 0)
        final_message_count = t;
    for (unsigned i = 0; i < final_message_count; ++i)
        high_partial[i] = 1;
    const unsigned high_message_shift =
        t + ((k - 1) / t) * t;
    run_case<Field>(
        ops, t, high_message_shift, true, bytes,
        high_partial, all_t, UINT64_C(0x4849474849464654), counts);

    std::vector<uint8_t> transmitted_high(t, 0);
    std::vector<uint8_t> requested_high(t, 0);
    for (unsigned i = 0; i < r; ++i)
    {
        transmitted_high[i] = 1;
        requested_high[i] = static_cast<uint8_t>((i % 3U) != 1U);
    }
    run_case<Field>(
        ops, t, 0, false, bytes, all_t, transmitted_high,
        UINT64_C(0x4849474850524546), counts);
    run_case<Field>(
        ops, t, 0, false, bytes, all_t, requested_high,
        UINT64_C(0x48494748484f4c45), counts);

    const unsigned p = ceil_power_of_two(k);
    const unsigned low_parent = ceil_power_of_two(p + r);
    require(low_parent <= Field::order(), "low profile exceeds field");
    std::vector<uint8_t> low_message(p, 0);
    std::vector<uint8_t> all_p(p, 1);
    for (unsigned i = 0; i < k; ++i)
        low_message[i] = 1;
    run_case<Field>(
        ops, p, 0, true, bytes, low_message, all_p,
        UINT64_C(0x4c4f574949464654), counts);

    unsigned final_recovery_count = r % p;
    if (final_recovery_count == 0)
        final_recovery_count = p;
    std::vector<uint8_t> low_final(p, 0);
    for (unsigned i = 0; i < final_recovery_count; ++i)
        low_final[i] = 1;
    const unsigned low_parity_shift = p + ((r - 1) / p) * p;
    run_case<Field>(
        ops, p, low_parity_shift, false, bytes, all_p, low_final,
        UINT64_C(0x4c4f575052454649), counts);

    std::vector<uint8_t> requested_low(p, 0);
    const unsigned first_block_count = r < p ? r : p;
    for (unsigned i = 0; i < first_block_count; ++i)
        requested_low[i] = static_cast<uint8_t>((i % 5U) == 0U || i + 1 == p);
    run_case<Field>(
        ops, p, p, false, bytes, all_p, requested_low,
        UINT64_C(0x4c4f57484f4c4559), counts);
}

bool same_plan(
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
        left.zero_outputs != right.zero_outputs ||
        left.operations.size() != right.operations.size() ||
        left.fused_four_starts != right.fused_four_starts ||
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
    return true;
}

uint16_t invalid_log_provider(const void*, uint32_t)
{
    return 256;
}

void test_invalid_plan_construction()
{
    leopard2_internal::PrunedTransformPlan plan;
    plan.size = 77;
    plan.shift = 33;
    plan.zero_multiplier_log = 19;
    plan.first_layer_multiplier_log = 18;
    plan.inverse = true;
    plan.input_mask.push_back(9);
    plan.output_mask.push_back(8);
    leopard2_internal::PrunedTransformOperation operation = { 4, 5, 6, 7 };
    plan.operations.push_back(operation);
    plan.fused_four_starts.push_back(1);
    plan.zero_outputs.push_back(3);
    plan.full_butterfly_count = 11;
    plan.one_output_butterflies = 12;
    plan.input_zero_specializations = 13;
    plan.zero_multiplier_butterflies = 14;
    plan.one_multiplier_butterflies = 15;
    const leopard2_internal::PrunedTransformPlan snapshot(plan);

    uint8_t valid[8] = { 1, 1, 1, 1, 1, 1, 1, 1 };
    uint8_t malformed[8] = { 1, 1, 1, 2, 1, 1, 1, 1 };
    require(!leopard::ff8::PreparePrunedTransformPlan(
            3, 0, false, valid, valid, plan) && same_plan(plan, snapshot),
        "invalid size changed caller plan");
    require(!leopard::ff8::PreparePrunedTransformPlan(
            8, 1, false, valid, valid, plan) && same_plan(plan, snapshot),
        "invalid shift changed caller plan");
    require(!leopard::ff8::PreparePrunedTransformPlan(
            8, 0, false, malformed, valid, plan) && same_plan(plan, snapshot),
        "invalid mask changed caller plan");
    require(!leopard2_internal::CompilePrunedTransformPlan(
            256, 255, 8, 0, false, valid, valid,
            invalid_log_provider, NULL, plan) && same_plan(plan, snapshot),
        "invalid multiplier changed caller plan");
}

void test_metrics()
{
    const unsigned size = 16;
    std::vector<uint8_t> input(size, 0);
    std::vector<uint8_t> output(size, 0);
    for (unsigned i = 0; i < 9; ++i)
        input[i] = 1;
    for (unsigned i = 0; i < 7; ++i)
        output[i] = 1;
    leopard2_internal::PrunedTransformPlan plan;
    require(leopard::ff8::PreparePrunedTransformPlan(
            size, 0, false, input.data(), output.data(), plan),
        "metrics plan construction failed");
    require(plan.full_butterfly_count == 32,
        "incorrect full radix-2 butterfly count");
    require(!plan.operations.empty() &&
            plan.operations.size() < plan.full_butterfly_count,
        "illustrative ragged plan did not prune work");
    require(plan.one_output_butterflies != 0 ||
            plan.input_zero_specializations != 0,
        "illustrative ragged plan omitted boundary specialization");
}

template<class Field>
void test_fused_four_descriptors()
{
    const unsigned size = 16;
    const std::vector<uint8_t> all(size, 1);
    for (unsigned inverse = 0; inverse < 2; ++inverse)
    {
        leopard2_internal::PrunedTransformPlan plan;
        require(Field::prepare(
                size, 0, inverse != 0, all.data(), all.data(), plan),
            "complete fused plan construction failed");
        require(plan.fused_four_starts.size() ==
                (size / 4) * (exact_log2(size) / 2),
            "complete plan omitted four-way layer groups");
        for (size_t i = 0; i < plan.fused_four_starts.size(); ++i)
        {
            require(plan.fused_four_starts[i] + 3 < plan.operations.size(),
                "fused descriptor exceeds operation schedule");
            if (i != 0)
                require(plan.fused_four_starts[i - 1] + 4 <=
                        plan.fused_four_starts[i],
                    "fused descriptors overlap or are unsorted");
        }
    }
}

void test_max_parent_plan_footprint(TestCounts& counts)
{
    const unsigned size = 65536;
    const std::vector<uint8_t> all(size, 1);
    leopard2_internal::PrunedTransformPlan plan;
    require(GF16::prepare(
            size, 0, false, all.data(), all.data(), plan),
        "maximum GF16 plan construction failed");
    require(plan.operations.size() ==
            static_cast<size_t>(size / 2) * exact_log2(size),
        "maximum GF16 plan lost full-transform operations");
    require(plan.fused_four_starts.size() ==
            static_cast<size_t>(size / 4) * (exact_log2(size) / 2),
        "maximum GF16 plan omitted all-level fusion");
    const uint64_t bytes = plan_storage_bytes(plan);
    require(bytes < UINT64_C(16) * 1024 * 1024,
        "maximum GF16 plan exceeds bounded metadata budget");
    if (bytes > counts.max_plan_bytes)
        counts.max_plan_bytes = bytes;
}

void test_shared_plan_concurrency(const leopard::backend::Ops& ops)
{
    const unsigned size = 64;
    const uint64_t bytes = 129;
    const std::vector<std::vector<uint8_t> > masks = make_masks(size);
    leopard2_internal::PrunedTransformPlan plan;
    require(leopard::ff8::PreparePrunedTransformPlan(
            size, 128, false, masks[3].data(), masks[2].data(), plan),
        "concurrent plan construction failed");

    const std::vector<std::vector<uint8_t> > initial = make_input(
        size, bytes, masks[3], UINT64_C(0x434f4e4355525245));
    std::vector<std::vector<uint8_t> > expected(initial);
    std::vector<void*> expected_pointers = pointers(expected);
    GF8::full(ops, false, bytes, size, 128, expected_pointers.data());

    std::atomic<bool> failed(false);
    std::vector<std::thread> threads;
    for (unsigned thread = 0; thread < 8; ++thread)
    {
        threads.push_back(std::thread([&, thread]() {
            for (unsigned iteration = 0; iteration < 16; ++iteration)
            {
                std::vector<std::vector<uint8_t> > actual(initial);
                std::vector<void*> actual_pointers = pointers(actual);
                if (!leopard2_internal::ExecutePrunedTransformPlan(
                        ops, bytes, plan, actual_pointers.data()))
                {
                    failed.store(true, std::memory_order_relaxed);
                    return;
                }
                for (unsigned coordinate = 0; coordinate < size; ++coordinate)
                    if (masks[2][coordinate] &&
                        actual[coordinate] != expected[coordinate])
                    {
                        (void)thread;
                        failed.store(true, std::memory_order_relaxed);
                        return;
                    }
            }
        }));
    }
    for (size_t i = 0; i < threads.size(); ++i)
        threads[i].join();
    require(!failed.load(std::memory_order_relaxed),
        "shared immutable plan failed concurrent execution");
}

void test_shared_source_plan_concurrency(const leopard::backend::Ops& ops)
{
    const unsigned size = 64;
    const uint64_t bytes = 129;
    const std::vector<uint8_t> input(size, 1);
    const std::vector<std::vector<uint8_t> > masks = make_masks(size);
    const std::vector<uint8_t>& output = masks[3];
    leopard2_internal::PrunedTransformPlan plan;
    require(leopard::ff8::PreparePrunedTransformPlan(
            size, 128, false, input.data(), output.data(), plan),
        "concurrent immutable-source plan construction failed");

    std::vector<std::vector<uint8_t> > source = make_input(
        size, bytes, input, UINT64_C(0x534f55524345504c));
    const std::vector<std::vector<uint8_t> > source_snapshot(source);
    std::vector<void*> source_values = source_pointers(source);
    std::vector<std::vector<uint8_t> > expected(source);
    std::vector<void*> expected_pointers = pointers(expected);
    GF8::full(ops, false, bytes, size, 128, expected_pointers.data());

    std::atomic<bool> failed(false);
    std::vector<std::thread> threads;
    for (unsigned thread = 0; thread < 8; ++thread)
    {
        threads.push_back(std::thread([&, thread]() {
            for (unsigned iteration = 0; iteration < 16; ++iteration)
            {
                std::vector<std::vector<uint8_t> > actual(
                    size, std::vector<uint8_t>(bytes, 0));
                std::vector<void*> actual_pointers = pointers(actual);
                if (!leopard2_internal::
                        ExecutePrunedForwardTransformPlanFromSources(
                            ops, bytes, plan, source_values.data(),
                            actual_pointers.data()))
                {
                    failed.store(true, std::memory_order_relaxed);
                    return;
                }
                for (unsigned coordinate = 0; coordinate < size; ++coordinate)
                    if (output[coordinate] &&
                        actual[coordinate] != expected[coordinate])
                    {
                        (void)thread;
                        failed.store(true, std::memory_order_relaxed);
                        return;
                    }
            }
        }));
    }
    for (size_t i = 0; i < threads.size(); ++i)
        threads[i].join();
    require(!failed.load(std::memory_order_relaxed),
        "shared immutable-source plan failed concurrent execution");
    require(source == source_snapshot,
        "concurrent immutable-source execution changed coefficients");
}

} // namespace

int main()
{
    try
    {
        require(leo_init() == Leopard_Success, "Leopard initialization failed");
        test_invalid_plan_construction();
        test_metrics();
        test_fused_four_descriptors<GF8>();
        test_fused_four_descriptors<GF16>();

        const leo2_backend requested[] = {
            LEO2_BACKEND_SCALAR,
            LEO2_BACKEND_SSSE3,
            LEO2_BACKEND_AVX2
        };
        const unsigned sizes8[] = { 2, 4, 8, 16, 64, 128, 256 };
        const unsigned sizes16[] = { 2, 4, 16, 64, 256, 1024 };
        const uint64_t bytes8[] = { 1, 7, 63, 64, 65, 129 };
        const uint64_t bytes16[] = { 2, 18, 62, 64, 66, 130 };
        TestCounts counts;
        test_max_parent_plan_footprint(counts);
        for (size_t i = 0; i < sizeof(requested) / sizeof(requested[0]); ++i)
        {
            leopard::backend::QualificationStatus status =
                leopard::backend::QualificationAvailable;
            const leopard::backend::Ops* ops =
                leopard::backend::GetQualifiedOps(requested[i], &status);
            if (!ops)
            {
                require(requested[i] != LEO2_BACKEND_SCALAR &&
                        status == leopard::backend::QualificationUnavailable,
                    "backend qualification failed unexpectedly");
                continue;
            }
            test_matrix<GF8>(
                *ops, sizes8, sizeof(sizes8) / sizeof(sizes8[0]),
                bytes8, sizeof(bytes8) / sizeof(bytes8[0]), counts);
            test_matrix<GF16>(
                *ops, sizes16, sizeof(sizes16) / sizeof(sizes16[0]),
                bytes16, sizeof(bytes16) / sizeof(bytes16[0]), counts);
            test_exhaustive_small_masks<GF8>(*ops, 17, counts);
            test_exhaustive_small_masks<GF16>(*ops, 18, counts);
            test_direct_oracle_masks<GF8>(*ops, counts);
            test_direct_oracle_masks<GF16>(*ops, counts);
            test_profile_masks<GF8>(*ops, 100, 30, 65, counts);
            test_profile_masks<GF16>(*ops, 1000, 200, 130, counts);
            test_profile_masks<GF8>(*ops, 17, 100, 129, counts);
            test_profile_masks<GF16>(*ops, 257, 700, 66, counts);
            test_shared_plan_concurrency(*ops);
            test_shared_source_plan_concurrency(*ops);
            ++counts.backends;
        }
        require(counts.backends != 0, "no backend was available");

        std::cout << "PASS pruned_transform"
                  << " backends=" << counts.backends
                  << " plans=" << counts.plans
                  << " executions=" << counts.executions
                  << " compared_bytes=" << counts.compared_bytes
                  << " direct_symbols=" << counts.direct_symbols
                  << " fused_four=" << counts.fused_four_descriptors
                  << " execution_steps=" << counts.execution_steps
                  << " max_plan_bytes=" << counts.max_plan_bytes
                  << std::endl;
        return 0;
    }
    catch (const std::exception& error)
    {
        std::cerr << "FAIL pruned_transform: " << error.what() << std::endl;
        return 1;
    }
}
