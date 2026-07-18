/*
    Copyright (c) 2017 Christopher A. Taylor.  All rights reserved.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the conditions in the Leopard-RS
    BSD license are met.  See the repository LICENSE file.
*/

#include "direct_oracle.h"
#include "Leopard2Direct.h"
#include "LeopardFF8.h"
#include "LeopardFF16.h"
#include "leopard2.h"

#include <algorithm>
#include <stdint.h>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <vector>

namespace {

using leopard2_test::BinaryField;
using leopard2_test::Element;
using leopard2_test::Matrix;
using leopard2_test::ProfileLayout;

typedef std::vector<std::vector<uint8_t> > Shards;

struct ProfileCase
{
    unsigned k;
    unsigned r;
    leo2_profile profile;
    leo2_field field;
};

struct Counts
{
    uint64_t profiles;
    uint64_t executions;
    uint64_t expected_symbols;
    uint64_t compared_bytes;
    uint64_t sparse_outputs;
    uint64_t generator_cross_checks;

    Counts()
        : profiles(0), executions(0), expected_symbols(0), compared_bytes(0)
        , sparse_outputs(0), generator_cross_checks(0)
    {
    }
};

class AlignedBuffer
{
public:
    explicit AlignedBuffer(size_t bytes)
        : pointer_(NULL), bytes_(bytes)
    {
        if (bytes_ != 0 && posix_memalign(
                &pointer_, leo2_scratch_alignment(), bytes_) != 0)
            pointer_ = NULL;
    }

    ~AlignedBuffer() { free(pointer_); }
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

leopard2_test::ProfileKind profile_kind(leo2_profile profile)
{
    return profile == LEO2_PROFILE_LEGACY_HIGH_V1
        ? leopard2_test::kLegacyHigh : leopard2_test::kLow;
}

Matrix barycentric_generator(
    const BinaryField& field,
    const ProfileLayout& layout)
{
    require(layout.parent_size <= field.order(),
        "profile exceeds oracle field");
    Matrix generator(layout.original_count + layout.recovery_count,
        std::vector<Element>(layout.original_count, 0));
    for (unsigned original = 0; original < layout.original_count; ++original)
        generator[original][original] = 1;

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

    for (unsigned parity = 0; parity < layout.recovery_count; ++parity)
    {
        const Element point = static_cast<Element>(
            layout.parity_coordinates[parity]);
        Element vanishing = 1;
        for (unsigned systematic = 0;
             systematic < layout.parent_dimension;
             ++systematic)
        {
            vanishing = field.multiply(vanishing, field.add(point,
                static_cast<Element>(layout.systematic_coordinates[systematic])));
        }
        for (unsigned original = 0; original < layout.original_count; ++original)
        {
            const Element difference = field.add(point,
                static_cast<Element>(layout.systematic_coordinates[original]));
            generator[layout.original_count + parity][original] =
                field.multiply(field.multiply(vanishing,
                    field.inverse(difference)), weights[original]);
        }
    }
    return generator;
}

void cross_check_generator(
    const BinaryField& field,
    leopard2_test::ProfileKind kind,
    unsigned k,
    unsigned r,
    Counts* counts)
{
    const ProfileLayout layout = leopard2_test::make_profile_layout(kind, k, r);
    require(barycentric_generator(field, layout) ==
            leopard2_test::direct_systematic_generator(field, layout),
        "barycentric and direct Lagrange generators differ");
    ++counts->generator_cross_checks;
}

size_t symbol_count(leo2_field field, size_t bytes)
{
    return field == LEO2_FIELD_GF8 ? bytes : bytes / 2;
}

Element read_symbol(
    const std::vector<uint8_t>& shard,
    leo2_field field,
    size_t symbol)
{
    if (field == LEO2_FIELD_GF8)
        return shard[symbol];
    const size_t total_symbols = shard.size() / 2;
    const size_t tile_symbol = symbol & 31u;
    const size_t tile_first = symbol - tile_symbol;
    const size_t tile_symbols = std::min<size_t>(32, total_symbols - tile_first);
    const size_t tile_byte = tile_first * 2;
    return static_cast<Element>(shard[tile_byte + tile_symbol] |
        (static_cast<unsigned>(shard[tile_byte + tile_symbols + tile_symbol]) << 8));
}

void write_symbol(
    std::vector<uint8_t>* shard,
    leo2_field field,
    size_t symbol,
    Element value)
{
    if (field == LEO2_FIELD_GF8)
    {
        (*shard)[symbol] = static_cast<uint8_t>(value);
        return;
    }
    const size_t total_symbols = shard->size() / 2;
    const size_t tile_symbol = symbol & 31u;
    const size_t tile_first = symbol - tile_symbol;
    const size_t tile_symbols = std::min<size_t>(32, total_symbols - tile_first);
    const size_t tile_byte = tile_first * 2;
    (*shard)[tile_byte + tile_symbol] = static_cast<uint8_t>(value);
    (*shard)[tile_byte + tile_symbols + tile_symbol] =
        static_cast<uint8_t>(value >> 8);
}

bool sparse_requested(unsigned parity, unsigned recovery_count)
{
    return parity == 0 || parity + 1 == recovery_count ||
        (parity % 7u) == 3u;
}

void compare_outputs(
    const Shards& actual,
    const Matrix& generator,
    const Shards& originals,
    const BinaryField& oracle_field,
    leo2_field field,
    bool sparse,
    Counts* counts)
{
    const size_t symbols = symbol_count(field, originals[0].size());
    const unsigned k = static_cast<unsigned>(originals.size());
    const unsigned r = static_cast<unsigned>(actual.size());
    std::vector<Element> message(k, 0);
    for (size_t symbol = 0; symbol < symbols; ++symbol)
    {
        for (unsigned original = 0; original < k; ++original)
            message[original] = read_symbol(originals[original], field, symbol);
        const std::vector<Element> expected =
            leopard2_test::matrix_vector_multiply(
                oracle_field, generator, message);
        for (unsigned parity = 0; parity < r; ++parity)
        {
            if (sparse && !sparse_requested(parity, r))
                continue;
            require(read_symbol(actual[parity], field, symbol) == expected[k + parity],
                "transform encoder differs from direct generator");
            ++counts->expected_symbols;
        }
    }
    for (unsigned parity = 0; parity < r; ++parity)
        if (!sparse || sparse_requested(parity, r))
        {
            counts->compared_bytes += actual[parity].size();
            if (sparse)
                ++counts->sparse_outputs;
        }
}

void run_profile(
    leo2_context* context,
    const ProfileCase& test,
    size_t bytes,
    Counts* counts)
{
    const BinaryField field = test.field == LEO2_FIELD_GF8
        ? leopard2_test::make_legacy_gf8()
        : leopard2_test::make_legacy_gf16();
    const ProfileLayout layout = leopard2_test::make_profile_layout(
        profile_kind(test.profile), test.k, test.r);
    const Matrix generator = barycentric_generator(field, layout);

    leo2_codec* codec = NULL;
    require_result(leo2_codec_create(context, test.k, test.r, test.profile,
        test.field, NULL, &codec), "transform codec create");
    require_result(leo2_test_codec_set_encode_mode(
        codec, LEO2_TEST_ENCODE_FORCE_TRANSFORM),
        "force transform encoder");
    require(leo2_codec_parent_count(codec) == layout.parent_size &&
            leo2_codec_padded_side(codec) == layout.padded_side &&
            leo2_codec_field(codec) == test.field,
        "transform profile introspection differs from direct layout");

    Shards originals(test.k, std::vector<uint8_t>(bytes, 0));
    uint64_t random = UINT64_C(0x4c324c434854524e) ^
        (static_cast<uint64_t>(test.k) << 32) ^
        (static_cast<uint64_t>(test.r) << 8) ^ bytes ^
        static_cast<uint64_t>(test.profile);
    const size_t symbols = symbol_count(test.field, bytes);
    for (unsigned original = 0; original < test.k; ++original)
        for (size_t symbol = 0; symbol < symbols; ++symbol)
            write_symbol(&originals[original], test.field, symbol,
                static_cast<Element>(splitmix64(&random) & (field.order() - 1)));
    std::vector<const void*> input(test.k, NULL);
    for (unsigned original = 0; original < test.k; ++original)
        input[original] = &originals[original][0];

    size_t scratch_bytes = 0;
    require_result(leo2_encode_scratch_size(codec, bytes, &scratch_bytes),
        "transform scratch query");
    AlignedBuffer scratch(scratch_bytes);
    require(scratch.valid(), "transform scratch allocation");

    for (unsigned sparse_value = 0; sparse_value < 2; ++sparse_value)
    {
        const bool sparse = sparse_value != 0;
        Shards parity(test.r, std::vector<uint8_t>(bytes, 0xa5));
        std::vector<void*> output(test.r, NULL);
        for (unsigned recovery = 0; recovery < test.r; ++recovery)
            if (!sparse || sparse_requested(recovery, test.r))
                output[recovery] = &parity[recovery][0];
        require_result(leo2_encode(codec, bytes, &input[0], &output[0],
            scratch.data(), scratch.size()), "transform encode");
        compare_outputs(parity, generator, originals, field, test.field,
            sparse, counts);
        ++counts->executions;
    }
    leo2_codec_destroy(codec);
    ++counts->profiles;
}

} // namespace

int main()
{
    try
    {
        Counts counts;
        cross_check_generator(leopard2_test::make_gf4(),
            leopard2_test::kLegacyHigh, 3, 2, &counts);
        cross_check_generator(leopard2_test::make_gf4(),
            leopard2_test::kLow, 3, 5, &counts);
        cross_check_generator(leopard2_test::make_legacy_gf8(),
            leopard2_test::kLegacyHigh, 5, 3, &counts);
        cross_check_generator(leopard2_test::make_legacy_gf8(),
            leopard2_test::kLow, 5, 11, &counts);

        leo2_context* context = NULL;
        require_result(leo2_context_create(NULL, &context),
            "transform context create");
        leopard::ff8::TestOnlyResetSparseEncodeCounts();
        leopard::ff16::TestOnlyResetSparseEncodeCounts();
        const ProfileCase cases[] = {
            { 3, 2, LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8 },
            { 5, 3, LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8 },
            { 9, 7, LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8 },
            { 17, 15, LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8 },
            { 31, 33, LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8 },
            { 65, 63, LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8 },
            { 129, 2, LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8 },
            { 191, 64, LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8 },
            { 240, 16, LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8 },
            { 248, 8, LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8 },
            { 2, 5, LEO2_PROFILE_LOW_V1, LEO2_FIELD_GF8 },
            { 3, 5, LEO2_PROFILE_LOW_V1, LEO2_FIELD_GF8 },
            { 5, 11, LEO2_PROFILE_LOW_V1, LEO2_FIELD_GF8 },
            { 7, 9, LEO2_PROFILE_LOW_V1, LEO2_FIELD_GF8 },
            { 9, 17, LEO2_PROFILE_LOW_V1, LEO2_FIELD_GF8 },
            { 15, 33, LEO2_PROFILE_LOW_V1, LEO2_FIELD_GF8 },
            { 17, 65, LEO2_PROFILE_LOW_V1, LEO2_FIELD_GF8 },
            { 31, 129, LEO2_PROFILE_LOW_V1, LEO2_FIELD_GF8 },
            { 32, 224, LEO2_PROFILE_LOW_V1, LEO2_FIELD_GF8 },
            { 64, 192, LEO2_PROFILE_LOW_V1, LEO2_FIELD_GF8 },
            { 257, 33, LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF16 },
            { 127, 129, LEO2_PROFILE_LOW_V1, LEO2_FIELD_GF16 }
        };
        for (size_t i = 0; i < sizeof(cases) / sizeof(cases[0]); ++i)
        {
            if (cases[i].field == LEO2_FIELD_GF8)
            {
                run_profile(context, cases[i], 1, &counts);
                run_profile(context, cases[i], 17, &counts);
            }
            else
            {
                run_profile(context, cases[i], 2, &counts);
                run_profile(context, cases[i], 66, &counts);
            }
        }
        const leopard::ff8::TestOnlySparseEncodeCounts sparse8 =
            leopard::ff8::TestOnlyGetSparseEncodeCounts();
        const leopard::ff16::TestOnlySparseEncodeCounts sparse16 =
            leopard::ff16::TestOnlyGetSparseEncodeCounts();
        require(sparse8.exact_blocks != 0 &&
                sparse8.retained_butterflies < sparse8.prefix_butterflies,
            "GF8 forced sparse path did not prune transform work");
        require(sparse16.exact_blocks != 0 &&
                sparse16.retained_butterflies < sparse16.prefix_butterflies,
            "GF16 forced sparse path did not prune transform work");
        leo2_context_destroy(context);
        std::cout << "leopard2 transform differential passed: profiles="
                  << counts.profiles << " executions=" << counts.executions
                  << " expected_symbols=" << counts.expected_symbols
                  << " compared_bytes=" << counts.compared_bytes
                  << " sparse_outputs=" << counts.sparse_outputs
                  << " sparse8_exact_blocks=" << sparse8.exact_blocks
                  << " sparse8_retained=" << sparse8.retained_butterflies
                  << " sparse8_prefix=" << sparse8.prefix_butterflies
                  << " sparse16_exact_blocks=" << sparse16.exact_blocks
                  << " sparse16_retained=" << sparse16.retained_butterflies
                  << " sparse16_prefix=" << sparse16.prefix_butterflies
                  << " generator_cross_checks=" << counts.generator_cross_checks
                  << std::endl;
        return 0;
    }
    catch (const std::exception& error)
    {
        std::cerr << "leopard2 transform differential failed: "
                  << error.what() << std::endl;
        return 1;
    }
}
