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
    ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
    LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
    CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
    SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
    INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
    CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
    ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
    POSSIBILITY OF SUCH DAMAGE.
*/

/*
    This compatibility matrix deliberately keeps the API/framing paths separate:
    Leopard2 receives compact application buffers, while the old leo_encode API
    receives explicitly packed 64-byte GF16 ALTMAP tiles.  The old API is the
    byte-for-byte compatibility oracle for wire layout, tails, and parity order;
    the direct-oracle tests separately provide arithmetic independence.
*/

#include "leopard.h"
#include "leopard2.h"

#include <algorithm>
#include <cstdint>
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

typedef std::vector<std::vector<uint8_t> > Shards;

struct Profile
{
    unsigned original_count;
    unsigned recovery_count;
    unsigned parent_count;
};

struct Counts
{
    uint64_t cases;
    uint64_t parity_bytes;
    uint64_t subset_parity_bytes;

    Counts()
        : cases(0)
        , parity_bytes(0)
        , subset_parity_bytes(0)
    {}
};

void require(bool condition, const std::string& message)
{
    if (!condition)
        throw std::runtime_error(message);
}

void require_result(leo2_result result, const std::string& operation)
{
    if (result != LEO2_SUCCESS)
    {
        std::ostringstream stream;
        stream << operation << ": " << leo2_result_string(result)
               << " (" << static_cast<int>(result) << ")";
        throw std::runtime_error(stream.str());
    }
}

struct AlignedBuffer
{
    void* data;
    size_t bytes;

    explicit AlignedBuffer(size_t size)
        : data(NULL)
        , bytes(size)
    {
        if (size != 0)
        {
#if defined(_MSC_VER)
            data = _aligned_malloc(size, leo2_scratch_alignment());
            if (!data)
                throw std::bad_alloc();
#else
            if (posix_memalign(&data, leo2_scratch_alignment(), size) != 0)
                throw std::bad_alloc();
#endif
            memset(data, 0, size);
        }
    }

    ~AlignedBuffer()
    {
#if defined(_MSC_VER)
        _aligned_free(data);
#else
        free(data);
#endif
    }

private:
    AlignedBuffer(const AlignedBuffer&);
    AlignedBuffer& operator=(const AlignedBuffer&);
};

uint64_t mix_seed(
    uint64_t seed,
    unsigned original_count,
    unsigned recovery_count,
    size_t bytes)
{
    seed ^= static_cast<uint64_t>(original_count) * UINT64_C(0x9e3779b97f4a7c15);
    seed ^= static_cast<uint64_t>(recovery_count) * UINT64_C(0xd1b54a32d192ed03);
    seed ^= static_cast<uint64_t>(bytes) * UINT64_C(0x94d049bb133111eb);
    return seed ? seed : UINT64_C(0xa0761d6478bd642f);
}

Shards make_originals(unsigned count, size_t bytes, uint64_t seed)
{
    Shards shards(count, std::vector<uint8_t>(bytes));
    uint64_t state = seed;
    for (unsigned shard = 0; shard < count; ++shard)
    {
        for (size_t offset = 0; offset < bytes; ++offset)
        {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            shards[shard][offset] = static_cast<uint8_t>(
                state + shard * 29u + offset * 131u);
        }
    }
    return shards;
}

std::vector<const void*> const_pointers(const Shards& shards)
{
    std::vector<const void*> pointers(shards.size());
    for (size_t i = 0; i < shards.size(); ++i)
        pointers[i] = &shards[i][0];
    return pointers;
}

std::vector<void*> mutable_pointers(Shards& shards)
{
    std::vector<void*> pointers(shards.size());
    for (size_t i = 0; i < shards.size(); ++i)
        pointers[i] = &shards[i][0];
    return pointers;
}

std::vector<uint8_t> compact_pack_gf16(const std::vector<uint8_t>& input)
{
    require(!input.empty() && (input.size() & 1u) == 0,
        "GF16 compact pack requires a positive even byte count");
    const size_t rounded = (input.size() + 63u) & ~static_cast<size_t>(63u);
    std::vector<uint8_t> output(rounded, 0);
    const size_t complete = input.size() & ~static_cast<size_t>(63u);
    if (complete != 0)
        memcpy(&output[0], &input[0], complete);
    const size_t symbols = (input.size() - complete) / 2;
    if (symbols != 0)
    {
        memcpy(&output[complete], &input[complete], symbols);
        memcpy(&output[complete + 32], &input[complete + symbols], symbols);
    }
    return output;
}

std::vector<uint8_t> compact_gather_gf16(
    const std::vector<uint8_t>& input,
    size_t bytes)
{
    const size_t rounded = (bytes + 63u) & ~static_cast<size_t>(63u);
    require(bytes != 0 && (bytes & 1u) == 0 && input.size() == rounded,
        "GF16 compact gather requires a positive even byte count");
    std::vector<uint8_t> output(bytes, 0);
    const size_t complete = bytes & ~static_cast<size_t>(63u);
    if (complete != 0)
        memcpy(&output[0], &input[0], complete);
    const size_t symbols = (bytes - complete) / 2;
    if (symbols != 0)
    {
        memcpy(&output[complete], &input[complete], symbols);
        memcpy(&output[complete + symbols], &input[complete + 32], symbols);
    }
    return output;
}

Shards encode_legacy_gf16(
    const Shards& original,
    unsigned recovery_count,
    size_t bytes)
{
    const size_t rounded = (bytes + 63u) & ~static_cast<size_t>(63u);
    Shards packed(original.size());
    for (size_t i = 0; i < original.size(); ++i)
        packed[i] = compact_pack_gf16(original[i]);
    std::vector<const void*> original_ptrs = const_pointers(packed);

    const unsigned work_count = leo_encode_work_count(
        static_cast<unsigned>(original.size()), recovery_count);
    require(work_count != 0, "legacy work-count query rejected a valid profile");
    Shards work(work_count, std::vector<uint8_t>(rounded, 0));
    std::vector<void*> work_ptrs = mutable_pointers(work);
    const LeopardResult result = leo_encode(
        rounded,
        static_cast<unsigned>(original.size()),
        recovery_count,
        work_count,
        &original_ptrs[0],
        &work_ptrs[0]);
    require(result == Leopard_Success,
        std::string("legacy GF16 encode failed: ") + leo_result_string(result));

    Shards recovery(recovery_count);
    for (unsigned i = 0; i < recovery_count; ++i)
        recovery[i] = compact_gather_gf16(work[i], bytes);
    return recovery;
}

void compare_parity(
    const std::vector<uint8_t>& actual,
    const std::vector<uint8_t>& expected,
    const Profile& profile,
    size_t bytes,
    size_t seed_index,
    size_t parity,
    const char* gate)
{
    require(actual.size() == expected.size(), "parity byte-count mismatch");
    for (size_t offset = 0; offset < actual.size(); ++offset)
    {
        if (actual[offset] != expected[offset])
        {
            std::ostringstream stream;
            stream << gate << " mismatch: K=" << profile.original_count
                   << " R=" << profile.recovery_count
                   << " N=" << profile.parent_count
                   << " bytes=" << bytes
                   << " seed=" << seed_index
                   << " parity=" << parity
                   << " offset=" << offset
                   << " actual=" << static_cast<unsigned>(actual[offset])
                   << " expected=" << static_cast<unsigned>(expected[offset]);
            throw std::runtime_error(stream.str());
        }
    }
}

void compare_shards(
    const Shards& actual,
    const Shards& expected,
    const Profile& profile,
    size_t bytes,
    size_t seed_index,
    const char* gate)
{
    require(actual.size() == expected.size(), "parity shard-count mismatch");
    for (size_t parity = 0; parity < actual.size(); ++parity)
    {
        compare_parity(actual[parity], expected[parity], profile, bytes,
            seed_index, parity, gate);
    }
}

std::vector<unsigned> subset_indices(unsigned recovery_count)
{
    std::vector<unsigned> indices;
    indices.push_back(0);
    indices.push_back(recovery_count / 2);
    indices.push_back(recovery_count - 1);
    std::sort(indices.begin(), indices.end());
    indices.erase(std::unique(indices.begin(), indices.end()), indices.end());
    return indices;
}

void run_case(
    leo2_context* context,
    const Profile& profile,
    size_t bytes,
    uint64_t seed,
    size_t seed_index,
    Counts* counts)
{
    leo2_codec* codec = NULL;
    require_result(leo2_codec_create(
        context,
        profile.original_count,
        profile.recovery_count,
        LEO2_PROFILE_LEGACY_HIGH_V1,
        LEO2_FIELD_GF16,
        NULL,
        &codec), "GF16 legacy-high codec create");
    require(codec != NULL, "codec create returned null");
    require(leo2_codec_profile(codec) == LEO2_PROFILE_LEGACY_HIGH_V1,
        "codec selected the wrong profile");
    require(leo2_codec_field(codec) == LEO2_FIELD_GF16,
        "codec selected the wrong field");
    require(leo2_codec_parent_count(codec) == profile.parent_count,
        "codec parent differs from the matrix specification");

    const Shards original = make_originals(
        profile.original_count,
        bytes,
        mix_seed(seed, profile.original_count, profile.recovery_count, bytes));
    const Shards expected = encode_legacy_gf16(
        original, profile.recovery_count, bytes);

    std::vector<const void*> original_ptrs = const_pointers(original);
    Shards actual(profile.recovery_count, std::vector<uint8_t>(bytes, 0));
    std::vector<void*> actual_ptrs = mutable_pointers(actual);
    size_t scratch_bytes = 0;
    require_result(leo2_encode_scratch_size(codec, bytes, &scratch_bytes),
        "GF16 legacy-high scratch query");
    AlignedBuffer scratch(scratch_bytes);
    require_result(leo2_encode(
        codec,
        bytes,
        &original_ptrs[0],
        &actual_ptrs[0],
        scratch.data,
        scratch.bytes), "GF16 legacy-high full encode");
    compare_shards(actual, expected, profile, bytes, seed_index, "full parity");

    const std::vector<unsigned> selected = subset_indices(profile.recovery_count);
    Shards subset(profile.recovery_count, std::vector<uint8_t>(bytes, 0xa5));
    std::vector<void*> subset_ptrs(profile.recovery_count, NULL);
    std::vector<uint8_t> requested(profile.recovery_count, 0);
    for (size_t i = 0; i < selected.size(); ++i)
    {
        subset_ptrs[selected[i]] = &subset[selected[i]][0];
        requested[selected[i]] = 1;
    }
    require_result(leo2_encode(
        codec,
        bytes,
        &original_ptrs[0],
        &subset_ptrs[0],
        scratch.data,
        scratch.bytes), "GF16 legacy-high subset encode");
    for (unsigned parity = 0; parity < profile.recovery_count; ++parity)
    {
        if (requested[parity])
        {
            compare_parity(subset[parity], expected[parity], profile, bytes,
                seed_index, parity, "subset parity");
            counts->subset_parity_bytes += bytes;
        }
        else
        {
            require(static_cast<size_t>(std::count(
                        subset[parity].begin(), subset[parity].end(),
                        static_cast<uint8_t>(0xa5))) == bytes,
                "an unrequested parity buffer was modified");
        }
    }

    ++counts->cases;
    counts->parity_bytes +=
        static_cast<uint64_t>(profile.recovery_count) * bytes;
    leo2_codec_destroy(codec);
}

} // namespace

int main()
{
    try
    {
        /*
            All entries have N > 256, so the old API independently selects its
            GF16 implementation.  Exact and just-over-boundary pairs exercise
            both sides of several parent-size transitions.
        */
        const Profile profiles[] = {
            { 129, 128, 512 },
            { 193, 64, 512 },
            { 225, 32, 512 },
            { 241, 16, 512 },
            { 249, 8, 512 },
            { 253, 4, 512 },
            { 255, 2, 512 },
            { 257, 33, 512 },
            { 384, 128, 512 },
            { 385, 128, 1024 },
            { 449, 64, 1024 },
            { 481, 32, 1024 },
            { 505, 8, 1024 },
            { 511, 2, 1024 },
            { 768, 256, 1024 },
            { 769, 256, 2048 },
            { 1000, 200, 2048 },
            { 1536, 512, 2048 },
            { 1537, 512, 4096 },
            { 4096, 512, 8192 }
        };
        const size_t byte_counts[] = { 2, 6, 34, 62, 66, 130 };
        const uint64_t seeds[] = {
            UINT64_C(0x0123456789abcdef),
            UINT64_C(0xd1b54a32d192ed03),
            UINT64_C(0x9e3779b97f4a7c15)
        };

        leo2_context_options options;
        memset(&options, 0, sizeof(options));
        options.struct_size = sizeof(options);
        options.backend = LEO2_BACKEND_AUTO;
        options.thread_count = 1;
        leo2_context* context = NULL;
        require_result(leo2_context_create(&options, &context), "context create");
        require(context != NULL, "context create returned null");

        Counts counts;
        for (size_t profile_i = 0;
             profile_i < sizeof(profiles) / sizeof(profiles[0]);
             ++profile_i)
        {
            for (size_t seed_i = 0; seed_i < sizeof(seeds) / sizeof(seeds[0]); ++seed_i)
            {
                for (size_t bytes_i = 0;
                     bytes_i < sizeof(byte_counts) / sizeof(byte_counts[0]);
                     ++bytes_i)
                {
                    run_case(context, profiles[profile_i], byte_counts[bytes_i],
                        seeds[seed_i], seed_i, &counts);
                }
            }
        }

        leo2_context_destroy(context);
        std::cout << "GF16 legacy-high encoder compatibility passed: profiles="
                  << sizeof(profiles) / sizeof(profiles[0])
                  << " cases=" << counts.cases
                  << " seeds=" << sizeof(seeds) / sizeof(seeds[0])
                  << " parity_bytes=" << counts.parity_bytes
                  << " subset_parity_bytes=" << counts.subset_parity_bytes
                  << std::endl;
        return 0;
    }
    catch (const std::exception& exception)
    {
        std::cerr << "GF16 legacy-high encoder matrix failed: "
                  << exception.what() << std::endl;
        return 1;
    }
}
