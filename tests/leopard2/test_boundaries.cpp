/*
    Copyright (c) 2017 Christopher A. Taylor.  All rights reserved.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the conditions in the Leopard-RS
    BSD license are met.  See the repository LICENSE file.
*/

#include "leopard2.h"

#include <algorithm>
#include <stdint.h>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <vector>

#if defined(_MSC_VER)
#include <malloc.h>
#endif

namespace {

typedef std::vector<std::vector<uint8_t> > Shards;

struct BoundaryCase
{
    uint32_t k;
    uint32_t r;
    leo2_profile profile;
};

struct Counts
{
    uint64_t cases;
    uint64_t executions;
    uint64_t recovered_shards;
    uint64_t compared_bytes;
    uint64_t gf8_cases;
    uint64_t gf16_cases;

    Counts()
        : cases(0), executions(0), recovered_shards(0), compared_bytes(0)
        , gf8_cases(0), gf16_cases(0)
    {
    }
};

class AlignedBuffer
{
public:
    explicit AlignedBuffer(size_t bytes)
        : pointer_(NULL)
        , bytes_(bytes)
    {
        if (bytes_ != 0)
        {
#if defined(_MSC_VER)
            pointer_ = _aligned_malloc(bytes_, leo2_scratch_alignment());
#else
            if (posix_memalign(
                    &pointer_, leo2_scratch_alignment(), bytes_) != 0)
                pointer_ = NULL;
#endif
        }
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

Shards make_shards(uint32_t count, size_t bytes, uint64_t seed)
{
    Shards shards(count, std::vector<uint8_t>(bytes));
    for (uint32_t shard = 0; shard < count; ++shard)
        for (size_t i = 0; i < bytes; ++i)
            shards[shard][i] = static_cast<uint8_t>(splitmix64(&seed) >> 56);
    return shards;
}

std::vector<const void*> const_pointers(const Shards& shards)
{
    std::vector<const void*> result(shards.size(), NULL);
    for (size_t i = 0; i < shards.size(); ++i)
        result[i] = shards[i].empty() ? NULL : &shards[i][0];
    return result;
}

std::vector<void*> mutable_pointers(Shards* shards)
{
    std::vector<void*> result(shards->size(), NULL);
    for (size_t i = 0; i < shards->size(); ++i)
        result[i] = (*shards)[i].empty() ? NULL : &(*shards)[i][0];
    return result;
}

leo2_codec* make_codec(
    leo2_context* context,
    const BoundaryCase& test,
    uint32_t flags)
{
    leo2_codec_options options;
    std::memset(&options, 0, sizeof(options));
    options.struct_size = sizeof(options);
    options.flags = flags;
    leo2_codec* codec = NULL;
    require_result(leo2_codec_create(context, test.k, test.r, test.profile,
        LEO2_FIELD_AUTO, &options, &codec), "boundary codec create");
    require(codec != NULL, "boundary codec create returned null");
    return codec;
}

Shards encode(
    const leo2_codec* codec,
    const Shards& originals,
    uint32_t recovery_count,
    size_t bytes)
{
    Shards recovery(recovery_count, std::vector<uint8_t>(bytes, 0));
    std::vector<const void*> input = const_pointers(originals);
    std::vector<void*> output = mutable_pointers(&recovery);
    size_t scratch_bytes = 0;
    require_result(leo2_encode_scratch_size(codec, bytes, &scratch_bytes),
        "boundary encode scratch");
    AlignedBuffer scratch(scratch_bytes);
    require(scratch.valid(), "boundary encode scratch allocation");
    require_result(leo2_encode(codec, bytes, &input[0], &output[0],
        scratch.data(), scratch.size()), "boundary encode");
    return recovery;
}

std::vector<uint32_t> missing_originals(uint32_t k, uint32_t losses)
{
    std::vector<uint32_t> result;
    result.reserve(losses);
    for (uint32_t i = 0; i < losses; ++i)
    {
        const uint32_t index = losses == 1
            ? k / 2
            : static_cast<uint32_t>(
                (static_cast<uint64_t>(i) * (k - 1)) / (losses - 1));
        if (result.empty() || result.back() != index)
            result.push_back(index);
    }
    require(result.size() == losses, "boundary loss selection duplicated an index");
    return result;
}

void execute_plan(
    leo2_decode_plan* plan,
    const Shards& originals,
    const Shards& parity,
    const std::vector<uint8_t>& original_present,
    const std::vector<uint8_t>& parity_present,
    size_t bytes,
    Shards* restored,
    Counts* counts)
{
    std::vector<const void*> original_input(originals.size(), NULL);
    std::vector<const void*> parity_input(parity.size(), NULL);
    std::vector<void*> restored_output(originals.size(), NULL);
    for (size_t i = 0; i < originals.size(); ++i)
    {
        if (original_present[i])
            original_input[i] = &originals[i][0];
        else
            restored_output[i] = &(*restored)[i][0];
    }
    for (size_t i = 0; i < parity.size(); ++i)
        if (parity_present[i])
            parity_input[i] = &parity[i][0];

    size_t scratch_bytes = 0;
    require_result(leo2_decode_plan_scratch_size(plan, bytes, &scratch_bytes),
        "boundary decode scratch");
    AlignedBuffer scratch(scratch_bytes);
    require(scratch.valid(), "boundary decode scratch allocation");
    for (unsigned repeat = 0; repeat < 2; ++repeat)
    {
        for (size_t i = 0; i < restored->size(); ++i)
            if (!original_present[i])
                std::fill((*restored)[i].begin(), (*restored)[i].end(),
                    static_cast<uint8_t>(0xa5u ^ repeat));
        require_result(leo2_decode_plan_execute(plan, bytes,
            &original_input[0], &parity_input[0], &restored_output[0],
            scratch.data(), scratch.size()), "boundary decode execute");
        for (size_t i = 0; i < restored->size(); ++i)
            if (!original_present[i])
            {
                require((*restored)[i] == originals[i],
                    "boundary recovered original mismatch");
                ++counts->recovered_shards;
                counts->compared_bytes += bytes;
            }
        ++counts->executions;
    }
}

void run_boundary_case(
    leo2_context* context,
    const BoundaryCase& test,
    size_t bytes,
    Counts* counts)
{
    leo2_codec* specialized = make_codec(context, test, 0);
    leo2_codec* generic = make_codec(
        context, test, LEO2_CODEC_FORCE_GENERIC_DECODE);
    require(leo2_codec_parent_count(specialized) ==
            leo2_codec_parent_count(generic) &&
        leo2_codec_padded_side(specialized) == leo2_codec_padded_side(generic) &&
        leo2_codec_field(specialized) == leo2_codec_field(generic),
        "generic flag changed boundary code identity");

    const leo2_field field = leo2_codec_field(specialized);
    require(field == LEO2_FIELD_GF8 || field == LEO2_FIELD_GF16,
        "boundary AUTO field was not resolved");
    if (field == LEO2_FIELD_GF8)
        ++counts->gf8_cases;
    else
        ++counts->gf16_cases;
    require(field == LEO2_FIELD_GF8 || (bytes & 1u) == 0,
        "boundary test selected an odd GF16 length");

    const uint64_t seed = UINT64_C(0x4c32424f554e4400) ^
        (static_cast<uint64_t>(test.k) << 32) ^
        (static_cast<uint64_t>(test.r) << 8) ^ bytes ^
        static_cast<uint64_t>(test.profile);
    const Shards originals = make_shards(test.k, bytes, seed);
    const Shards parity = encode(specialized, originals, test.r, bytes);

    const uint32_t losses = std::min(std::min(test.k, test.r), 4u);
    const std::vector<uint32_t> missing = missing_originals(test.k, losses);
    std::vector<uint8_t> original_present(test.k, 1);
    std::vector<uint8_t> parity_present(test.r, 0);
    for (size_t i = 0; i < missing.size(); ++i)
        original_present[missing[i]] = 0;
    // Leave exactly L parity equations.  Together with K-L surviving
    // originals this is the maximum-erasure boundary of the public code.
    for (uint32_t i = 0; i < losses; ++i)
        parity_present[i] = 1;

    leo2_decode_plan* specialized_plan = NULL;
    leo2_decode_plan* generic_plan = NULL;
    require_result(leo2_decode_plan_create(specialized,
        &original_present[0], &parity_present[0], &specialized_plan),
        "boundary specialized plan");
    require_result(leo2_decode_plan_create(generic,
        &original_present[0], &parity_present[0], &generic_plan),
        "boundary generic plan");

    Shards specialized_restored(test.k, std::vector<uint8_t>(bytes, 0));
    Shards generic_restored(test.k, std::vector<uint8_t>(bytes, 0));
    execute_plan(specialized_plan, originals, parity, original_present,
        parity_present, bytes, &specialized_restored, counts);
    execute_plan(generic_plan, originals, parity, original_present,
        parity_present, bytes, &generic_restored, counts);
    for (size_t i = 0; i < missing.size(); ++i)
        require(specialized_restored[missing[i]] ==
                generic_restored[missing[i]],
            "specialized/generic boundary differential mismatch");

    leo2_decode_plan_destroy(generic_plan);
    leo2_decode_plan_destroy(specialized_plan);
    leo2_codec_destroy(generic);
    leo2_codec_destroy(specialized);
    ++counts->cases;
}

void test_field_limits(leo2_context* context)
{
    leo2_codec* codec = NULL;
    require(leo2_codec_create(context, 225, 30,
        LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8, NULL, &codec) ==
        LEO2_INVALID_COUNTS, "GF8 accepted an inflated parent over 256");
    require(codec == NULL, "failed GF8 field-limit create returned a codec");
    require_result(leo2_codec_create(context, 225, 30,
        LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF16, NULL, &codec),
        "GF16 inflated-parent create");
    require(leo2_codec_parent_count(codec) == 512,
        "GF16 inflated-parent count mismatch");
    leo2_codec_destroy(codec);
}

} // namespace

int main()
{
    try
    {
        leo2_context* context = NULL;
        require_result(leo2_context_create(NULL, &context),
            "boundary context create");
        test_field_limits(context);

        const BoundaryCase cases[] = {
            { 1, 1, LEO2_PROFILE_LEGACY_HIGH_V1 },
            { 1, 8, LEO2_PROFILE_LOW_V1 },
            { 3, 5, LEO2_PROFILE_LOW_V1 },
            { 7, 9, LEO2_PROFILE_LOW_V1 },
            { 8, 248, LEO2_PROFILE_LOW_V1 },
            { 9, 7, LEO2_PROFILE_LEGACY_HIGH_V1 },
            { 15, 17, LEO2_PROFILE_LOW_V1 },
            { 16, 240, LEO2_PROFILE_LOW_V1 },
            { 17, 15, LEO2_PROFILE_LEGACY_HIGH_V1 },
            { 31, 33, LEO2_PROFILE_LOW_V1 },
            { 32, 224, LEO2_PROFILE_LOW_V1 },
            { 33, 31, LEO2_PROFILE_LEGACY_HIGH_V1 },
            { 63, 65, LEO2_PROFILE_LOW_V1 },
            { 64, 192, LEO2_PROFILE_LOW_V1 },
            { 65, 63, LEO2_PROFILE_LEGACY_HIGH_V1 },
            { 100, 156, LEO2_PROFILE_LOW_V1 },
            { 100, 156, LEO2_PROFILE_LEGACY_HIGH_V1 },
            { 127, 129, LEO2_PROFILE_LOW_V1 },
            { 128, 128, LEO2_PROFILE_LEGACY_HIGH_V1 },
            { 129, 1, LEO2_PROFILE_LEGACY_HIGH_V1 },
            { 129, 100, LEO2_PROFILE_LEGACY_HIGH_V1 },
            { 191, 64, LEO2_PROFILE_LEGACY_HIGH_V1 },
            { 192, 64, LEO2_PROFILE_LEGACY_HIGH_V1 },
            { 223, 32, LEO2_PROFILE_LEGACY_HIGH_V1 },
            { 224, 32, LEO2_PROFILE_LEGACY_HIGH_V1 },
            { 225, 30, LEO2_PROFILE_LEGACY_HIGH_V1 },
            { 239, 16, LEO2_PROFILE_LEGACY_HIGH_V1 },
            { 240, 16, LEO2_PROFILE_LEGACY_HIGH_V1 },
            { 241, 15, LEO2_PROFILE_LEGACY_HIGH_V1 },
            { 247, 8, LEO2_PROFILE_LEGACY_HIGH_V1 },
            { 248, 8, LEO2_PROFILE_LEGACY_HIGH_V1 },
            { 249, 7, LEO2_PROFILE_LEGACY_HIGH_V1 },
            { 255, 1, LEO2_PROFILE_LEGACY_HIGH_V1 },
            { 256, 256, LEO2_PROFILE_LOW_V1 },
            { 256, 256, LEO2_PROFILE_LEGACY_HIGH_V1 },
            { 257, 255, LEO2_PROFILE_LEGACY_HIGH_V1 },
            { 511, 513, LEO2_PROFILE_LOW_V1 },
            { 512, 512, LEO2_PROFILE_LEGACY_HIGH_V1 },
            { 513, 511, LEO2_PROFILE_LEGACY_HIGH_V1 },
            { 1000, 200, LEO2_PROFILE_LEGACY_HIGH_V1 },
            { 4096, 512, LEO2_PROFILE_LEGACY_HIGH_V1 }
        };

        Counts counts;
        for (size_t case_i = 0; case_i < sizeof(cases) / sizeof(cases[0]);
             ++case_i)
        {
            leo2_codec* probe = make_codec(context, cases[case_i], 0);
            const leo2_field field = leo2_codec_field(probe);
            leo2_codec_destroy(probe);
            const size_t gf8_bytes[] = { 1, 17, 65 };
            const size_t gf16_bytes[] = { 2, 66 };
            if (field == LEO2_FIELD_GF8)
                for (size_t i = 0; i < sizeof(gf8_bytes) / sizeof(gf8_bytes[0]); ++i)
                    run_boundary_case(context, cases[case_i], gf8_bytes[i], &counts);
            else
                for (size_t i = 0; i < sizeof(gf16_bytes) / sizeof(gf16_bytes[0]); ++i)
                    run_boundary_case(context, cases[case_i], gf16_bytes[i], &counts);
        }

        leo2_context_destroy(context);
        std::cout << "leopard2 boundary differential passed: cases="
                  << counts.cases << " executions=" << counts.executions
                  << " recovered_shards=" << counts.recovered_shards
                  << " compared_bytes=" << counts.compared_bytes
                  << " gf8_cases=" << counts.gf8_cases
                  << " gf16_cases=" << counts.gf16_cases << std::endl;
        return 0;
    }
    catch (const std::exception& error)
    {
        std::cerr << "leopard2 boundary differential failed: "
                  << error.what() << std::endl;
        return 1;
    }
}
