/*
    Copyright (c) 2017 Christopher A. Taylor.  All rights reserved.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the conditions in the Leopard-RS
    BSD license are met.  See the repository LICENSE file.
*/

#include "direct_oracle.h"
#include "leopard2.h"

#include <algorithm>
#include <stdint.h>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <new>
#include <sstream>
#include <stdexcept>
#include <vector>

#if defined(_MSC_VER)
#include <malloc.h>
#endif

namespace {

using leopard2_test::BinaryField;
using leopard2_test::Element;
using leopard2_test::Matrix;
using leopard2_test::ProfileLayout;

struct LowCase
{
    unsigned k;
    unsigned r;
    size_t bytes;
};

struct Counts
{
    uint64_t profiles;
    uint64_t executions;
    uint64_t partial_masks;
    uint64_t requested_outputs;
    uint64_t direct_rows;
    uint64_t direct_coefficients;
    uint64_t compared_symbols;
    uint64_t compared_bytes;
    uint64_t dense_crosscheck_rows;
    uint64_t dense_crosscheck_coefficients;
    unsigned largest_parent;

    Counts()
        : profiles(0), executions(0), partial_masks(0), requested_outputs(0)
        , direct_rows(0), direct_coefficients(0), compared_symbols(0)
        , compared_bytes(0), dense_crosscheck_rows(0)
        , dense_crosscheck_coefficients(0), largest_parent(0)
    {
    }
};

class AlignedBuffer
{
public:
    explicit AlignedBuffer(size_t bytes)
        : pointer_(NULL), bytes_(bytes)
    {
        if (bytes_ == 0)
            return;
#if defined(_MSC_VER)
        pointer_ = _aligned_malloc(bytes_, leo2_scratch_alignment());
#else
        if (posix_memalign(&pointer_, leo2_scratch_alignment(), bytes_) != 0)
            pointer_ = NULL;
#endif
    }

    ~AlignedBuffer()
    {
#if defined(_MSC_VER)
        _aligned_free(pointer_);
#else
        free(pointer_);
#endif
    }
    bool valid() const { return bytes_ == 0 || pointer_ != NULL; }
    void* data() { return pointer_; }
    size_t size() const { return bytes_; }

private:
    AlignedBuffer(const AlignedBuffer&);
    AlignedBuffer& operator=(const AlignedBuffer&);
    void* pointer_;
    size_t bytes_;
};

void require(bool condition, const std::string& message)
{
    if (!condition)
        throw std::runtime_error(message);
}

void require_result(leo2_result result, const char* operation)
{
    if (result != LEO2_SUCCESS)
    {
        std::ostringstream stream;
        stream << operation << ": " << leo2_result_string(result);
        throw std::runtime_error(stream.str());
    }
}

uint64_t splitmix64(uint64_t* state)
{
    uint64_t value = (*state += UINT64_C(0x9e3779b97f4a7c15));
    value = (value ^ (value >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
    value = (value ^ (value >> 27)) * UINT64_C(0x94d049bb133111eb);
    return value ^ (value >> 31);
}

Element read_gf16_symbol(const std::vector<uint8_t>& shard, size_t symbol)
{
    const size_t total_symbols = shard.size() / 2;
    const size_t tile_symbol = symbol & 31u;
    const size_t tile_first = symbol - tile_symbol;
    const size_t tile_symbols = std::min<size_t>(32, total_symbols - tile_first);
    const size_t tile_byte = tile_first * 2;
    return static_cast<Element>(shard[tile_byte + tile_symbol] |
        (static_cast<unsigned>(
            shard[tile_byte + tile_symbols + tile_symbol]) << 8));
}

void write_gf16_symbol(
    std::vector<uint8_t>* shard,
    size_t symbol,
    Element value)
{
    const size_t total_symbols = shard->size() / 2;
    const size_t tile_symbol = symbol & 31u;
    const size_t tile_first = symbol - tile_symbol;
    const size_t tile_symbols = std::min<size_t>(32, total_symbols - tile_first);
    const size_t tile_byte = tile_first * 2;
    (*shard)[tile_byte + tile_symbol] = static_cast<uint8_t>(value);
    (*shard)[tile_byte + tile_symbols + tile_symbol] =
        static_cast<uint8_t>(value >> 8);
}

std::vector<Element> low_systematic_weights(
    const BinaryField& field,
    const ProfileLayout& layout)
{
    require(layout.kind == leopard2_test::kLow,
        "direct-row oracle requires the low profile");
    std::vector<Element> weights(layout.original_count, 0);
    for (unsigned original = 0; original < layout.original_count; ++original)
    {
        const Element x = static_cast<Element>(
            layout.systematic_coordinates[original]);
        Element denominator = 1;
        for (unsigned other = 0; other < layout.parent_dimension; ++other)
        {
            const Element y = static_cast<Element>(
                layout.systematic_coordinates[other]);
            if (x != y)
                denominator = field.multiply(denominator, field.add(x, y));
        }
        weights[original] = field.inverse(denominator);
    }
    return weights;
}

std::vector<Element> low_direct_parity_row(
    const BinaryField& field,
    const ProfileLayout& layout,
    const std::vector<Element>& weights,
    unsigned parity_index)
{
    require(weights.size() == layout.original_count,
        "direct-row weight count mismatch");
    require(parity_index < layout.recovery_count,
        "direct-row parity index out of range");
    const Element point = static_cast<Element>(
        layout.parity_coordinates[parity_index]);
    Element vanishing = 1;
    for (unsigned systematic = 0;
         systematic < layout.parent_dimension;
         ++systematic)
    {
        vanishing = field.multiply(vanishing, field.add(point,
            static_cast<Element>(layout.systematic_coordinates[systematic])));
    }

    std::vector<Element> row(layout.original_count, 0);
    for (unsigned original = 0; original < layout.original_count; ++original)
    {
        const Element difference = field.add(point,
            static_cast<Element>(layout.systematic_coordinates[original]));
        row[original] = field.multiply(field.multiply(vanishing,
            field.inverse(difference)), weights[original]);
    }
    return row;
}

std::vector<unsigned> requested_rows(unsigned r)
{
    std::vector<unsigned> requested;
    requested.push_back(0);
    requested.push_back(r / 4);
    requested.push_back(r / 2);
    requested.push_back(r - 1);
    std::sort(requested.begin(), requested.end());
    requested.erase(std::unique(requested.begin(), requested.end()),
        requested.end());
    return requested;
}

std::vector<std::vector<unsigned> > partial_masks(unsigned r)
{
    std::vector<std::vector<unsigned> > masks;
    masks.push_back(std::vector<unsigned>(1, 0));
    masks.push_back(std::vector<unsigned>(1, r - 1));
    masks.push_back(requested_rows(r));
    return masks;
}

size_t selected_row_index(
    const std::vector<unsigned>& selected,
    unsigned parity)
{
    const std::vector<unsigned>::const_iterator found =
        std::lower_bound(selected.begin(), selected.end(), parity);
    require(found != selected.end() && *found == parity,
        "partial mask referenced an unbuilt oracle row");
    return static_cast<size_t>(found - selected.begin());
}

void cross_check_dense_generator(
    const BinaryField& field,
    unsigned k,
    unsigned r,
    Counts* counts)
{
    const ProfileLayout layout = leopard2_test::make_profile_layout(
        leopard2_test::kLow, k, r);
    const Matrix generator = leopard2_test::direct_systematic_generator(
        field, layout);
    const std::vector<Element> weights = low_systematic_weights(field, layout);
    for (unsigned parity = 0; parity < r; ++parity)
    {
        const std::vector<Element> row = low_direct_parity_row(
            field, layout, weights, parity);
        require(row == generator[k + parity],
            "selected direct row differs from direct generator");
        ++counts->dense_crosscheck_rows;
        counts->dense_crosscheck_coefficients += row.size();
    }
}

void run_case(
    leo2_context* context,
    const BinaryField& field,
    const LowCase& test,
    Counts* counts)
{
    require(test.r > test.k, "GF16 low boundary case is not low rate");
    require(test.bytes != 0 && (test.bytes & 1u) == 0,
        "GF16 low boundary byte count is not even");
    const ProfileLayout layout = leopard2_test::make_profile_layout(
        leopard2_test::kLow, test.k, test.r);
    require(layout.parent_size <= field.order(),
        "GF16 low boundary profile exceeds the field");

    leo2_codec* codec = NULL;
    require_result(leo2_codec_create(context, test.k, test.r,
        LEO2_PROFILE_LOW_V1, LEO2_FIELD_GF16, NULL, &codec),
        "GF16 low boundary codec create");
    require(leo2_codec_parent_count(codec) == layout.parent_size &&
            leo2_codec_padded_side(codec) == layout.padded_side,
        "GF16 low boundary introspection mismatch");

    std::vector<std::vector<uint8_t> > originals(
        test.k, std::vector<uint8_t>(test.bytes, 0));
    uint64_t random = UINT64_C(0x4c4f574746313652) ^
        (static_cast<uint64_t>(test.k) << 32) ^
        (static_cast<uint64_t>(test.r) << 8) ^ test.bytes;
    const size_t symbols = test.bytes / 2;
    for (unsigned original = 0; original < test.k; ++original)
        for (size_t symbol = 0; symbol < symbols; ++symbol)
            write_gf16_symbol(&originals[original], symbol,
                static_cast<Element>(splitmix64(&random)));
    std::vector<const void*> input(test.k, NULL);
    for (unsigned original = 0; original < test.k; ++original)
        input[original] = &originals[original][0];

    const std::vector<unsigned> selected = requested_rows(test.r);
    const std::vector<Element> weights = low_systematic_weights(field, layout);
    Matrix direct_rows;
    for (size_t i = 0; i < selected.size(); ++i)
    {
        direct_rows.push_back(low_direct_parity_row(
            field, layout, weights, selected[i]));
        ++counts->direct_rows;
        counts->direct_coefficients += test.k;
    }

    size_t scratch_bytes = 0;
    require_result(leo2_encode_scratch_size(codec, test.bytes, &scratch_bytes),
        "GF16 low boundary scratch query");
    AlignedBuffer scratch(scratch_bytes);
    require(scratch.valid(), "GF16 low boundary scratch allocation");

    const std::vector<std::vector<unsigned> > masks = partial_masks(test.r);
    for (size_t mask_index = 0; mask_index < masks.size(); ++mask_index)
    {
        const std::vector<unsigned>& mask = masks[mask_index];
        std::vector<std::vector<uint8_t> > actual(
            mask.size(), std::vector<uint8_t>(test.bytes, 0xa5));
        std::vector<void*> output(test.r, NULL);
        for (size_t output_index = 0; output_index < mask.size(); ++output_index)
            output[mask[output_index]] = &actual[output_index][0];
        require_result(leo2_encode(codec, test.bytes, &input[0], &output[0],
            scratch.data(), scratch.size()), "GF16 low boundary encode");

        for (size_t output_index = 0; output_index < mask.size(); ++output_index)
        {
            const std::vector<Element>& row = direct_rows[selected_row_index(
                selected, mask[output_index])];
            for (size_t symbol = 0; symbol < symbols; ++symbol)
            {
                Element expected = 0;
                for (unsigned original = 0; original < test.k; ++original)
                {
                    expected ^= field.multiply(row[original],
                        read_gf16_symbol(originals[original], symbol));
                }
                require(read_gf16_symbol(actual[output_index], symbol) == expected,
                    "GF16 low transform differs from selected direct row");
                ++counts->compared_symbols;
            }
            counts->compared_bytes += test.bytes;
            ++counts->requested_outputs;
        }
        ++counts->executions;
        ++counts->partial_masks;
    }

    counts->largest_parent = std::max(counts->largest_parent,
        layout.parent_size);
    ++counts->profiles;
    leo2_codec_destroy(codec);
}

} // namespace

int main()
{
    try
    {
        Counts counts;
        const BinaryField field = leopard2_test::make_legacy_gf16();
        cross_check_dense_generator(field, 3, 5, &counts);
        cross_check_dense_generator(field, 5, 11, &counts);
        cross_check_dense_generator(field, 7, 9, &counts);

        leo2_context_options options;
        memset(&options, 0, sizeof(options));
        options.struct_size = sizeof(options);
        options.backend = LEO2_BACKEND_AUTO;
        options.thread_count = 1;
        leo2_context* context = NULL;
        require_result(leo2_context_create(&options, &context),
            "GF16 low boundary context create");
        const LowCase cases[] = {
            { 3, 5, 2 },
            { 7, 9, 6 },
            { 15, 17, 66 },
            { 17, 33, 2 },
            { 31, 33, 6 },
            { 33, 65, 66 },
            { 63, 65, 2 },
            { 65, 129, 6 },
            { 127, 129, 66 },
            { 129, 257, 2 },
            { 255, 257, 6 },
            { 257, 513, 66 }
        };
        for (size_t i = 0; i < sizeof(cases) / sizeof(cases[0]); ++i)
            run_case(context, field, cases[i], &counts);
        leo2_context_destroy(context);

        std::cout << "leopard2 GF16 low direct-row differential passed: profiles="
                  << counts.profiles << " executions=" << counts.executions
                  << " partial_masks=" << counts.partial_masks
                  << " requested_outputs=" << counts.requested_outputs
                  << " direct_rows=" << counts.direct_rows
                  << " direct_coefficients=" << counts.direct_coefficients
                  << " compared_symbols=" << counts.compared_symbols
                  << " compared_bytes=" << counts.compared_bytes
                  << " dense_crosscheck_rows=" << counts.dense_crosscheck_rows
                  << " dense_crosscheck_coefficients="
                  << counts.dense_crosscheck_coefficients
                  << " largest_parent=" << counts.largest_parent << std::endl;
        return 0;
    }
    catch (const std::exception& error)
    {
        std::cerr << "leopard2 GF16 low direct-row differential failed: "
                  << error.what() << std::endl;
        return 1;
    }
}
