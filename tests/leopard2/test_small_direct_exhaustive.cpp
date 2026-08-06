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

#include "direct_oracle.h"
#include "Leopard2Dispatch.h"
#include "leopard2.h"

#include <cerrno>
#include <cstdlib>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <vector>

#ifndef LEO2_EXPERIMENT_GF8_SMALL_DIRECT_MODE
#define LEO2_EXPERIMENT_GF8_SMALL_DIRECT_MODE 0
#endif

namespace {

static const uint64_t kExpectedPatternCount = UINT64_C(1982812);
static const uint64_t kFnvOffset = UINT64_C(14695981039346656037);
static const uint64_t kFnvPrime = UINT64_C(1099511628211);

using leopard2_test::BinaryField;
using leopard2_test::Matrix;
using leopard2_test::ProfileLayout;

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
        stream << operation << ": " << leo2_result_string(result)
               << " (" << static_cast<int>(result) << ")";
        throw std::runtime_error(stream.str());
    }
}

uint32_t parse_u32(const char* text, const char* label)
{
    errno = 0;
    char* end = NULL;
    const unsigned long long value = std::strtoull(text, &end, 10);
    if (errno != 0 || !end || *end != '\0' || end == text ||
        value > std::numeric_limits<uint32_t>::max())
        throw std::invalid_argument(std::string("invalid ") + label);
    return static_cast<uint32_t>(value);
}

bool next_combination(std::vector<unsigned>& values, unsigned universe)
{
    if (values.empty())
        return false;
    for (size_t reverse = values.size(); reverse != 0; --reverse)
    {
        const size_t index = reverse - 1;
        const unsigned maximum = universe -
            static_cast<unsigned>(values.size() - index);
        if (values[index] == maximum)
            continue;
        ++values[index];
        for (size_t suffix = index + 1; suffix < values.size(); ++suffix)
            values[suffix] = values[suffix - 1] + 1;
        return true;
    }
    return false;
}

void initialize_combination(std::vector<unsigned>& values)
{
    for (size_t i = 0; i < values.size(); ++i)
        values[i] = static_cast<unsigned>(i);
}

uint64_t mix_u64(uint64_t digest, uint64_t value)
{
    for (unsigned byte = 0; byte < 8; ++byte)
    {
        digest ^= static_cast<uint8_t>(value >> (byte * 8));
        digest *= kFnvPrime;
    }
    return digest;
}

struct Counts
{
    uint64_t total_patterns;
    uint64_t assigned_patterns;
    uint64_t recovered_shards;
    uint64_t verified_basis_symbols;
    uint64_t digest;
    uint64_t ordinal_digest;

    Counts()
        : total_patterns(0)
        , assigned_patterns(0)
        , recovered_shards(0)
        , verified_basis_symbols(0)
        , digest(kFnvOffset)
        , ordinal_digest(kFnvOffset)
    {
    }
};

void verify_pattern(
    leo2_context* context,
    leo2_codec* codec,
    const Matrix& generator,
    unsigned k,
    unsigned r,
    const std::vector<unsigned>& missing,
    const std::vector<unsigned>& selected_parities,
    Counts& counts)
{
    std::vector<uint8_t> original_present(k, 1);
    std::vector<uint8_t> recovery_present(r, 0);
    for (size_t i = 0; i < missing.size(); ++i)
        original_present[missing[i]] = 0;
    for (size_t i = 0; i < selected_parities.size(); ++i)
        recovery_present[selected_parities[i]] = 1;

    leo2_decode_plan* plan = NULL;
    require_result(leo2_decode_plan_create(codec, original_present.data(),
        recovery_present.data(), &plan), "exhaustive plan create");
    try
    {
        leopard2_internal::DecodePathInfo path;
        require_result(leopard2_internal::GetDecodePlanPathInfo(
            plan, k, false, &path), "exhaustive path introspection");
        const bool production_full_loss_transform =
            LEO2_EXPERIMENT_GF8_SMALL_DIRECT_MODE == 0 &&
            missing.size() == k && k >= 7;
        if (production_full_loss_transform)
        {
            require(path.path != leopard2_internal::kDecodePathDirect &&
                    path.rule != leopard2_internal::kDecodeRuleDirect,
                "production full-loss transform control selected direct repair");
        }
        else
        {
            require(path.path == leopard2_internal::kDecodePathDirect &&
                    path.rule == leopard2_internal::kDecodeRuleDirect,
                "eligible loss pattern did not select direct repair");
        }
        const leopard2_internal::DirectRepairExecutor expected_executor =
            LEO2_EXPERIMENT_GF8_SMALL_DIRECT_MODE != 1
                ? leopard2_internal::kDirectRepairExecutorSourceMajor
                : leopard2_internal::kDirectRepairExecutorOutputMajor;
        if (!production_full_loss_transform)
        {
            require(path.direct_executor == expected_executor,
                "eligible loss pattern selected the wrong direct executor");
        }

        size_t scratch_bytes = 0;
        require_result(leo2_decode_plan_scratch_size(
            plan, k, &scratch_bytes), "exhaustive scratch query");
        std::vector<uint8_t> scratch(scratch_bytes + 63, 0);
        void* scratch_pointer = NULL;
        if (scratch_bytes != 0)
        {
            const uintptr_t raw = reinterpret_cast<uintptr_t>(scratch.data());
            scratch_pointer = reinterpret_cast<void*>((raw + 63u) & ~uintptr_t(63u));
        }

        /* Pack all K coordinate-basis messages into K byte lanes.  A single
           decode therefore proves the complete K-column linear recovery map
           for this missing/original and selected/parity coefficient matrix. */
        std::vector<std::vector<uint8_t> > originals(
            k, std::vector<uint8_t>(k, 0));
        std::vector<std::vector<uint8_t> > parities(
            r, std::vector<uint8_t>(k, 0));
        std::vector<std::vector<uint8_t> > restored(
            k, std::vector<uint8_t>(k, 0xa5));
        for (unsigned row = 0; row < k; ++row)
            originals[row][row] = 1;
        for (unsigned parity = 0; parity < r; ++parity)
            for (unsigned column = 0; column < k; ++column)
                parities[parity][column] = static_cast<uint8_t>(
                    generator[k + parity][column]);

        std::vector<const void*> original_input(k, NULL);
        std::vector<const void*> recovery_input(r, NULL);
        std::vector<void*> output(k, NULL);
        for (unsigned original = 0; original < k; ++original)
            if (original_present[original])
                original_input[original] = originals[original].data();
        for (unsigned parity = 0; parity < r; ++parity)
            if (recovery_present[parity])
                recovery_input[parity] = parities[parity].data();
        for (size_t i = 0; i < missing.size(); ++i)
            output[missing[i]] = restored[missing[i]].data();

        require_result(leo2_decode_plan_execute(plan, k,
            original_input.data(), recovery_input.data(), output.data(),
            scratch_pointer, scratch_bytes), "exhaustive direct execution");
        for (size_t i = 0; i < missing.size(); ++i)
        {
            const unsigned original = missing[i];
            for (unsigned column = 0; column < k; ++column)
            {
                const uint8_t expected = original == column ? 1 : 0;
                if (restored[original][column] != expected)
                {
                    std::ostringstream stream;
                    stream << "basis recovery mismatch K=" << k
                           << " R=" << r << " L=" << missing.size()
                           << " original=" << original
                           << " column=" << column;
                    throw std::runtime_error(stream.str());
                }
            }
        }

        ++counts.assigned_patterns;
        counts.recovered_shards += missing.size();
        counts.verified_basis_symbols += missing.size() * k;
        counts.digest = mix_u64(counts.digest, k);
        counts.digest = mix_u64(counts.digest, r);
        counts.digest = mix_u64(counts.digest, missing.size());
        for (size_t i = 0; i < missing.size(); ++i)
            counts.digest = mix_u64(counts.digest, missing[i]);
        for (size_t i = 0; i < selected_parities.size(); ++i)
            counts.digest = mix_u64(counts.digest, selected_parities[i]);
        counts.ordinal_digest = mix_u64(
            counts.ordinal_digest, counts.total_patterns);
    }
    catch (...)
    {
        leo2_decode_plan_destroy(plan);
        throw;
    }
    leo2_decode_plan_destroy(plan);
    (void)context;
}

void run_shard(uint32_t shard_index, uint32_t shard_count, Counts& counts)
{
    require(LEO2_EXPERIMENT_GF8_SMALL_DIRECT_MODE >= 0 &&
            LEO2_EXPERIMENT_GF8_SMALL_DIRECT_MODE <= 2,
        "exhaustive small-direct verifier requires mode 0, 1, or 2");
    require(shard_count != 0 && shard_index < shard_count,
        "shard index must be smaller than nonzero shard count");

    leo2_context_options context_options;
    std::memset(&context_options, 0, sizeof(context_options));
    context_options.struct_size = sizeof(context_options);
    context_options.backend = LEO2_BACKEND_AVX2;
    context_options.thread_count = 1;
    leo2_context* context = NULL;
    require_result(leo2_context_create(&context_options, &context),
        "explicit AVX2 context create");
    try
    {
        require(leo2_context_backend(context) == LEO2_BACKEND_AVX2,
            "explicit AVX2 context resolved a different backend");
        const BinaryField field = leopard2_test::make_legacy_gf8();
        for (unsigned k = 5; k <= 16; ++k)
            for (unsigned r = 5; r <= 8; ++r)
            {
                const ProfileLayout layout =
                    leopard2_test::make_profile_layout(
                        leopard2_test::kLegacyHigh, k, r);
                const Matrix generator =
                    leopard2_test::direct_systematic_generator(field, layout);
                leo2_codec* codec = NULL;
                leo2_codec_options codec_options;
                std::memset(&codec_options, 0, sizeof(codec_options));
                codec_options.struct_size = sizeof(codec_options);
                require_result(leo2_codec_create(context, k, r,
                    LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8,
                    &codec_options, &codec), "exhaustive codec create");
                try
                {
                    const unsigned maximum_loss = k < r ? k : r;
                    for (unsigned losses = 5;
                         losses <= maximum_loss; ++losses)
                    {
                        std::vector<unsigned> missing(losses);
                        initialize_combination(missing);
                        do
                        {
                            std::vector<unsigned> parities(losses);
                            initialize_combination(parities);
                            do
                            {
                                if (counts.total_patterns % shard_count ==
                                        shard_index)
                                {
                                    verify_pattern(context, codec, generator,
                                        k, r, missing, parities, counts);
                                }
                                ++counts.total_patterns;
                            }
                            while (next_combination(parities, r));
                        }
                        while (next_combination(missing, k));
                    }
                }
                catch (...)
                {
                    leo2_codec_destroy(codec);
                    throw;
                }
                leo2_codec_destroy(codec);
            }
        require(counts.total_patterns == kExpectedPatternCount,
            "exhaustive pattern count changed");
        const uint64_t expected_assigned =
            kExpectedPatternCount / shard_count +
            (shard_index < kExpectedPatternCount % shard_count ? 1 : 0);
        require(counts.assigned_patterns == expected_assigned,
            "deterministic shard assigned the wrong pattern count");
    }
    catch (...)
    {
        leo2_context_destroy(context);
        throw;
    }
    leo2_context_destroy(context);
}

} // namespace

int main(int argc, char** argv)
{
    try
    {
        uint32_t shard_index = 0;
        uint32_t shard_count = 1;
        for (int i = 1; i < argc; ++i)
        {
            require(i + 1 < argc, "missing command-line option value");
            if (std::strcmp(argv[i], "--shard-index") == 0)
                shard_index = parse_u32(argv[++i], "shard index");
            else if (std::strcmp(argv[i], "--shard-count") == 0)
                shard_count = parse_u32(argv[++i], "shard count");
            else
                throw std::invalid_argument(
                    std::string("unknown option: ") + argv[i]);
        }

        Counts counts;
        run_shard(shard_index, shard_count, counts);
        std::cout << "{\n"
                  << "  \"schema\": \"leopard2-small-direct-exhaustive/v1\",\n"
                  << "  \"mode\": "
                  << LEO2_EXPERIMENT_GF8_SMALL_DIRECT_MODE << ",\n"
                  << "  \"shard_index\": " << shard_index << ",\n"
                  << "  \"shard_count\": " << shard_count << ",\n"
                  << "  \"total_patterns\": " << counts.total_patterns << ",\n"
                  << "  \"assigned_patterns\": "
                  << counts.assigned_patterns << ",\n"
                  << "  \"recovered_shards\": "
                  << counts.recovered_shards << ",\n"
                  << "  \"verified_basis_symbols\": "
                  << counts.verified_basis_symbols << ",\n"
                  << "  \"basis_seed\": 0,\n"
                  << "  \"assignment\": \"global_ordinal_mod_shard_count\",\n"
                  << "  \"digest_fnv1a64\": \""
                  << std::hex << std::setw(16) << std::setfill('0')
                  << counts.digest << "\",\n"
                  << "  \"ordinal_digest_fnv1a64\": \""
                  << std::setw(16) << std::setfill('0')
                  << counts.ordinal_digest << "\"\n"
                  << "}\n";
        return 0;
    }
    catch (const std::exception& error)
    {
        std::cerr << "small-direct exhaustive verification failed: "
                  << error.what() << std::endl;
        return 1;
    }
}
