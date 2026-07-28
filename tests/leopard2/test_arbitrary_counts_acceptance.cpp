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
    Focused arbitrary-count acceptance for explicit LOW_V1 with K > R.  AUTO
    normally selects the legacy-high profile in this region, so these legal low
    profiles need an explicit production test.  The K=5,R=3 case exhausts all
    K-of-(K+R) public survivor subsets against the independent direct oracle.
*/

#include "direct_oracle.h"
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

using leopard2_test::BinaryField;
using leopard2_test::Element;
using leopard2_test::Matrix;
using leopard2_test::ProfileLayout;

typedef std::vector<uint8_t> Bytes;
typedef std::vector<Bytes> Shards;

struct Counts
{
    uint64_t profiles;
    uint64_t exhaustive_subsets;
    uint64_t specialized_executions;
    uint64_t generic_executions;
    uint64_t direct_parity_symbols;
    uint64_t direct_recovered_symbols;
    uint64_t restored_shards;
    uint64_t parity_rebuilds;
    uint64_t parity_bytes;
    uint64_t surplus_received;
    uint64_t rejected_counts;

    Counts()
        : profiles(0)
        , exhaustive_subsets(0)
        , specialized_executions(0)
        , generic_executions(0)
        , direct_parity_symbols(0)
        , direct_recovered_symbols(0)
        , restored_shards(0)
        , parity_rebuilds(0)
        , parity_bytes(0)
        , surplus_received(0)
        , rejected_counts(0)
    {}
};

void require(bool condition, const std::string& message)
{
    if (!condition)
        throw std::runtime_error(message);
}

void require_result(
    leo2_result actual,
    leo2_result expected,
    const std::string& operation)
{
    if (actual == expected)
        return;
    std::ostringstream stream;
    stream << operation << ": expected " << leo2_result_string(expected)
           << " (" << static_cast<int>(expected) << "), received "
           << leo2_result_string(actual) << " (" << static_cast<int>(actual) << ')';
    throw std::runtime_error(stream.str());
}

void require_success(leo2_result result, const std::string& operation)
{
    require_result(result, LEO2_SUCCESS, operation);
}

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
        if (!data_)
            throw std::bad_alloc();
#else
        if (posix_memalign(&data_, leo2_scratch_alignment(), bytes_) != 0)
            throw std::bad_alloc();
#endif
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
    size_t bytes() const { return bytes_; }

private:
    AlignedBuffer(const AlignedBuffer&);
    AlignedBuffer& operator=(const AlignedBuffer&);

    void* data_;
    size_t bytes_;
};

unsigned popcount(unsigned value)
{
    unsigned count = 0;
    while (value != 0)
    {
        count += value & 1u;
        value >>= 1;
    }
    return count;
}

uint64_t splitmix64(uint64_t* state)
{
    uint64_t value = (*state += UINT64_C(0x9e3779b97f4a7c15));
    value = (value ^ (value >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
    value = (value ^ (value >> 27)) * UINT64_C(0x94d049bb133111eb);
    return value ^ (value >> 31);
}

Shards make_originals(unsigned count, size_t bytes, uint64_t seed)
{
    Shards result(count, Bytes(bytes, 0));
    for (unsigned shard = 0; shard < count; ++shard)
    {
        for (size_t offset = 0; offset < bytes; ++offset)
        {
            result[shard][offset] = static_cast<uint8_t>(
                splitmix64(&seed) + shard * 29u + offset * 131u);
        }
    }
    return result;
}

std::vector<const void*> const_pointers(const Shards& shards)
{
    std::vector<const void*> result(shards.size(), NULL);
    for (size_t i = 0; i < shards.size(); ++i)
        result[i] = &shards[i][0];
    return result;
}

std::vector<void*> mutable_pointers(Shards* shards)
{
    std::vector<void*> result(shards->size(), NULL);
    for (size_t i = 0; i < shards->size(); ++i)
        result[i] = &(*shards)[i][0];
    return result;
}

leo2_codec* create_codec(
    leo2_context* context,
    unsigned k,
    unsigned r,
    leo2_field field,
    uint32_t flags)
{
    leo2_codec_options options;
    memset(&options, 0, sizeof(options));
    options.struct_size = sizeof(options);
    options.flags = flags;
    leo2_codec* codec = NULL;
    require_success(leo2_codec_create(context, k, r, LEO2_PROFILE_LOW_V1,
        field, &options, &codec), "low codec create");
    require(codec != NULL, "low codec create returned null");
    require(leo2_codec_profile(codec) == LEO2_PROFILE_LOW_V1,
        "explicit low codec selected a different profile");
    return codec;
}

Shards encode(
    const leo2_codec* codec,
    const Shards& originals,
    unsigned recovery_count,
    size_t bytes)
{
    Shards recovery(recovery_count, Bytes(bytes, 0));
    std::vector<const void*> input = const_pointers(originals);
    std::vector<void*> output = mutable_pointers(&recovery);
    size_t scratch_bytes = 0;
    require_success(leo2_encode_scratch_size(codec, bytes, &scratch_bytes),
        "encode scratch query");
    AlignedBuffer scratch(scratch_bytes);
    require_success(leo2_encode(codec, bytes, &input[0], &output[0],
        scratch.data(), scratch.bytes()), "encode");
    return recovery;
}

Shards execute_decode(
    const leo2_codec* codec,
    const Shards& originals,
    const Shards& recovery,
    const std::vector<uint8_t>& original_present,
    const std::vector<uint8_t>& recovery_present,
    size_t bytes)
{
    leo2_decode_plan* plan = NULL;
    require_success(leo2_decode_plan_create(codec, &original_present[0],
        &recovery_present[0], &plan), "decode plan create");

    std::vector<const void*> original_input = const_pointers(originals);
    std::vector<const void*> recovery_input = const_pointers(recovery);
    Shards restored(originals.size(), Bytes(bytes, 0xa5));
    std::vector<void*> output(originals.size(), NULL);
    for (size_t i = 0; i < originals.size(); ++i)
    {
        if (!original_present[i])
        {
            original_input[i] = NULL;
            output[i] = &restored[i][0];
        }
    }
    for (size_t i = 0; i < recovery.size(); ++i)
        if (!recovery_present[i])
            recovery_input[i] = NULL;

    size_t scratch_bytes = 0;
    require_success(leo2_decode_plan_scratch_size(plan, bytes, &scratch_bytes),
        "decode scratch query");
    AlignedBuffer scratch(scratch_bytes);
    require_success(leo2_decode_plan_execute(plan, bytes, &original_input[0],
        &recovery_input[0], &output[0], scratch.data(), scratch.bytes()),
        "decode execute");
    leo2_decode_plan_destroy(plan);
    return restored;
}

Shards combine_originals(
    const Shards& originals,
    const Shards& restored,
    const std::vector<uint8_t>& original_present)
{
    Shards complete = originals;
    for (size_t i = 0; i < originals.size(); ++i)
        if (!original_present[i])
            complete[i] = restored[i];
    return complete;
}

void require_recovered(
    const Shards& originals,
    const Shards& specialized,
    const Shards& generic,
    const std::vector<uint8_t>& original_present,
    Counts* counts)
{
    for (size_t i = 0; i < originals.size(); ++i)
    {
        if (!original_present[i])
        {
            require(specialized[i] == originals[i],
                "specialized low decoder restored the wrong original");
            require(generic[i] == originals[i],
                "generic low decoder restored the wrong original");
            require(generic[i] == specialized[i],
                "specialized and generic low recovery differ");
            ++counts->restored_shards;
        }
    }
}

void require_reencode(
    const leo2_codec* codec,
    const Shards& complete,
    const Shards& expected_recovery,
    size_t bytes,
    Counts* counts)
{
    const Shards rebuilt = encode(codec, complete,
        static_cast<unsigned>(expected_recovery.size()), bytes);
    require(rebuilt == expected_recovery,
        "parity after recovery differs from a fresh encode");
    ++counts->parity_rebuilds;
    counts->parity_bytes +=
        static_cast<uint64_t>(expected_recovery.size()) * bytes;
}

void test_exhaustive_low_k_greater_r(
    leo2_context* context,
    const BinaryField& field,
    Counts* counts)
{
    const unsigned k = 5;
    const unsigned r = 3;
    const size_t bytes = 17;
    leo2_codec* specialized = create_codec(context, k, r, LEO2_FIELD_GF8,
        LEO2_CODEC_FORCE_SPECIALIZED_DECODE);
    leo2_codec* generic = create_codec(context, k, r, LEO2_FIELD_GF8,
        LEO2_CODEC_FORCE_GENERIC_DECODE);
    require(leo2_codec_padded_side(specialized) == 8 &&
            leo2_codec_parent_count(specialized) == 16 &&
            leo2_codec_field(specialized) == LEO2_FIELD_GF8,
        "LOW(5,3) coordinate identity is wrong");

    const Shards originals = make_originals(
        k, bytes, UINT64_C(0x4c4f573500030011));
    const Shards recovery = encode(specialized, originals, r, bytes);
    require(encode(generic, originals, r, bytes) == recovery,
        "decoder selection changed LOW(5,3) parity");

    const ProfileLayout layout = leopard2_test::make_profile_layout(
        leopard2_test::kLow, k, r);
    const Matrix generator = leopard2_test::direct_systematic_generator(field, layout);
    std::vector<std::vector<Element> > codewords(bytes);
    for (size_t offset = 0; offset < bytes; ++offset)
    {
        std::vector<Element> message(k, 0);
        for (unsigned i = 0; i < k; ++i)
            message[i] = originals[i][offset];
        codewords[offset] = leopard2_test::matrix_vector_multiply(
            field, generator, message);
        for (unsigned i = 0; i < r; ++i)
        {
            require(recovery[i][offset] == codewords[offset][k + i],
                "LOW(5,3) parity differs from the direct generator");
            ++counts->direct_parity_symbols;
        }
    }

    const unsigned public_count = k + r;
    const unsigned mask_limit = 1u << public_count;
    for (unsigned mask = 0; mask < mask_limit; ++mask)
    {
        if (popcount(mask) != k)
            continue;
        std::vector<uint8_t> original_present(k, 0);
        std::vector<uint8_t> recovery_present(r, 0);
        std::vector<unsigned> received_rows;
        for (unsigned row = 0; row < public_count; ++row)
        {
            if ((mask & (1u << row)) == 0)
                continue;
            received_rows.push_back(row);
            if (row < k)
                original_present[row] = 1;
            else
                recovery_present[row - k] = 1;
        }

        const Shards specialized_restored = execute_decode(specialized,
            originals, recovery, original_present, recovery_present, bytes);
        const Shards generic_restored = execute_decode(generic,
            originals, recovery, original_present, recovery_present, bytes);
        ++counts->specialized_executions;
        ++counts->generic_executions;
        require_recovered(originals, specialized_restored, generic_restored,
            original_present, counts);

        for (size_t offset = 0; offset < bytes; ++offset)
        {
            std::vector<Element> received_values;
            for (size_t i = 0; i < received_rows.size(); ++i)
                received_values.push_back(codewords[offset][received_rows[i]]);
            const std::vector<Element> direct = leopard2_test::direct_recover(
                field, generator, received_rows, received_values);
            for (unsigned i = 0; i < k; ++i)
            {
                if (!original_present[i])
                {
                    require(specialized_restored[i][offset] == direct[i],
                        "LOW(5,3) recovery differs from direct interpolation");
                    ++counts->direct_recovered_symbols;
                }
            }
        }

        const Shards complete = combine_originals(
            originals, specialized_restored, original_present);
        require_reencode(specialized, complete, recovery, bytes, counts);
        ++counts->exhaustive_subsets;
    }
    require(counts->exhaustive_subsets == 56,
        "LOW(5,3) did not enumerate all 56 K-survivor subsets");
    ++counts->profiles;
    leo2_codec_destroy(generic);
    leo2_codec_destroy(specialized);
}

void run_transition_profile(
    leo2_context* context,
    unsigned k,
    size_t bytes,
    leo2_field expected_field,
    unsigned expected_padded,
    unsigned expected_parent,
    const std::vector<unsigned>& missing_originals,
    const std::vector<unsigned>& missing_recovery,
    Counts* counts)
{
    const unsigned r = 100;
    leo2_codec* specialized = create_codec(context, k, r, LEO2_FIELD_AUTO,
        LEO2_CODEC_FORCE_SPECIALIZED_DECODE);
    leo2_codec* generic = create_codec(context, k, r, LEO2_FIELD_AUTO,
        LEO2_CODEC_FORCE_GENERIC_DECODE);
    require(leo2_codec_padded_side(specialized) == expected_padded &&
            leo2_codec_parent_count(specialized) == expected_parent &&
            leo2_codec_field(specialized) == expected_field,
        "low field-transition coordinate identity is wrong");
    require(leo2_codec_padded_side(generic) == expected_padded &&
            leo2_codec_parent_count(generic) == expected_parent &&
            leo2_codec_field(generic) == expected_field,
        "generic flag changed low field-transition identity");

    const Shards originals = make_originals(k, bytes,
        UINT64_C(0x4c32415242434e54) + k * 257u);
    const Shards recovery = encode(specialized, originals, r, bytes);
    require(encode(generic, originals, r, bytes) == recovery,
        "decoder selection changed low field-transition parity");

    std::vector<uint8_t> original_present(k, 1);
    std::vector<uint8_t> recovery_present(r, 1);
    for (size_t i = 0; i < missing_originals.size(); ++i)
        original_present[missing_originals[i]] = 0;
    for (size_t i = 0; i < missing_recovery.size(); ++i)
        recovery_present[missing_recovery[i]] = 0;
    const unsigned received_count = k + r -
        static_cast<unsigned>(missing_originals.size() + missing_recovery.size());
    require(received_count > k, "transition case does not have surplus receives");
    counts->surplus_received += received_count - k;

    const Shards specialized_restored = execute_decode(specialized,
        originals, recovery, original_present, recovery_present, bytes);
    const Shards generic_restored = execute_decode(generic,
        originals, recovery, original_present, recovery_present, bytes);
    ++counts->specialized_executions;
    ++counts->generic_executions;
    require_recovered(originals, specialized_restored, generic_restored,
        original_present, counts);
    const Shards complete = combine_originals(
        originals, specialized_restored, original_present);
    require_reencode(specialized, complete, recovery, bytes, counts);

    ++counts->profiles;
    leo2_codec_destroy(generic);
    leo2_codec_destroy(specialized);
}

void test_field_transition(leo2_context* context, Counts* counts)
{
    const std::vector<unsigned> missing_originals_128 = { 0, 31, 64, 127 };
    const std::vector<unsigned> missing_originals_129 = { 0, 32, 64, 96, 128 };
    const std::vector<unsigned> missing_recovery = { 0, 1, 50, 99 };
    run_transition_profile(context, 128, 17, LEO2_FIELD_GF8, 128, 256,
        missing_originals_128, missing_recovery, counts);
    run_transition_profile(context, 129, 66, LEO2_FIELD_GF16, 256, 512,
        missing_originals_129, missing_recovery, counts);

    leo2_codec* rejected = reinterpret_cast<leo2_codec*>(
        static_cast<uintptr_t>(1));
    require_result(leo2_codec_create(context, 129, 100, LEO2_PROFILE_LOW_V1,
        LEO2_FIELD_GF8, NULL, &rejected), LEO2_INVALID_COUNTS,
        "LOW(129,100) explicit GF8 field inflation");
    require(rejected == NULL, "rejected GF8 codec did not clear its output");
    ++counts->rejected_counts;
}

void verify_expected_backend(const leo2_context* context)
{
    const char* expected = std::getenv("LEO2_EXPECT_BACKEND");
    if (!expected || expected[0] == '\0')
        return;
    leo2_backend backend = LEO2_BACKEND_AUTO;
    if (std::strcmp(expected, "scalar") == 0)
        backend = LEO2_BACKEND_SCALAR;
    else if (std::strcmp(expected, "ssse3") == 0)
        backend = LEO2_BACKEND_SSSE3;
    else if (std::strcmp(expected, "avx2") == 0)
        backend = LEO2_BACKEND_AVX2;
    else if (std::strcmp(expected, "avx512") == 0)
        backend = LEO2_BACKEND_AVX512;
    else if (std::strcmp(expected, "gfni") == 0 ||
             std::strcmp(expected, "avx2-gfni") == 0)
        backend = LEO2_BACKEND_GFNI;
    else
        throw std::runtime_error("invalid LEO2_EXPECT_BACKEND value");
    require(leo2_context_backend(context) == backend,
        "forced build selected the wrong runtime backend");
}

} // namespace

int main()
{
    try
    {
        leo2_context_options options;
        memset(&options, 0, sizeof(options));
        options.struct_size = sizeof(options);
        options.backend = LEO2_BACKEND_AUTO;
        options.thread_count = 1;
        leo2_context* context = NULL;
        require_success(leo2_context_create(&options, &context), "context create");
        require(context != NULL, "context create returned null");
        verify_expected_backend(context);

        Counts counts;
        const BinaryField gf8 = leopard2_test::make_legacy_gf8();
        test_exhaustive_low_k_greater_r(context, gf8, &counts);
        test_field_transition(context, &counts);
        leo2_context_destroy(context);

        std::cout << "arbitrary-count acceptance passed: profiles=" << counts.profiles
                  << " exhaustive_subsets=" << counts.exhaustive_subsets
                  << " specialized_executions=" << counts.specialized_executions
                  << " generic_executions=" << counts.generic_executions
                  << " direct_parity_symbols=" << counts.direct_parity_symbols
                  << " direct_recovered_symbols=" << counts.direct_recovered_symbols
                  << " restored_shards=" << counts.restored_shards
                  << " parity_rebuilds=" << counts.parity_rebuilds
                  << " parity_bytes=" << counts.parity_bytes
                  << " surplus_received=" << counts.surplus_received
                  << " rejected_counts=" << counts.rejected_counts
                  << std::endl;
        return 0;
    }
    catch (const std::exception& exception)
    {
        std::cerr << "arbitrary-count acceptance failed: "
                  << exception.what() << std::endl;
        return 1;
    }
}
