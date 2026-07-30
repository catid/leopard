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

#include "Leopard2Backend.h"
#include "Leopard2Direct.h"
#include "LeopardFF8.h"
#include "leopard.h"
#include "leopard2.h"

#include <algorithm>
#include <atomic>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <memory>
#include <new>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

#if defined(_MSC_VER)
#include <malloc.h>
#elif defined(__linux__)
#include <sys/stat.h>
#include <unistd.h>
#endif

#ifndef LEO2_ENABLE_TEST_HOOKS
#error "The AUTO encode-backend test requires LEO2_ENABLE_TEST_HOOKS"
#endif

namespace {

typedef std::vector<uint8_t> Bytes;
typedef std::vector<Bytes> Shards;

void require(bool condition, const char* message)
{
    if (!condition)
        throw std::runtime_error(message);
}

void require_result(
    leo2_result actual, leo2_result expected, const char* message)
{
    if (actual != expected)
        throw std::runtime_error(std::string(message) + ": " +
            leo2_result_string(actual));
}

void test_balanced_execution_tile_geometry()
{
    static const size_t requested_tile_bytes = 32U * 1024U;
    static const size_t requested_tile_count =
        requested_tile_bytes / leo2_scratch_alignment();
    static const size_t pass_counts[] = { 1, 2, 63, 64, 65, 129, 513 };
    for (size_t count_i = 0;
         count_i < sizeof(pass_counts) / sizeof(pass_counts[0]); ++count_i)
    {
        const size_t expected_count = pass_counts[count_i];
        for (size_t final_pass_tiles = 1;
             final_pass_tiles <= requested_tile_count;
             ++final_pass_tiles)
        {
            const size_t total_tiles =
                (expected_count - 1U) * requested_tile_count +
                final_pass_tiles;
            const size_t aligned_bytes =
                total_tiles * leo2_scratch_alignment();
            size_t execution_tile_count = 0;
            size_t maximum_pass_bytes = 0;
            require_result(leo2_test_balanced_execution_tiles(
                    aligned_bytes, requested_tile_bytes,
                    &execution_tile_count, &maximum_pass_bytes),
                LEO2_SUCCESS, "balanced execution-tile geometry");
            require(execution_tile_count == expected_count,
                "balanced execution-tile count differs from ceiling division");

            size_t remaining_tiles = total_tiles;
            size_t reference_maximum_tiles = 0;
            for (size_t pass = 0; pass < expected_count; ++pass)
            {
                const size_t passes_left = expected_count - pass;
                const size_t pass_tiles =
                    remaining_tiles / passes_left +
                    (remaining_tiles % passes_left != 0);
                reference_maximum_tiles = std::max(
                    reference_maximum_tiles, pass_tiles);
                remaining_tiles -= pass_tiles;
            }
            require(remaining_tiles == 0,
                "balanced execution reference did not consume every tile");
            require(maximum_pass_bytes ==
                    reference_maximum_tiles * leo2_scratch_alignment(),
                "balanced execution scratch is smaller than a distributed pass");
        }
    }

    size_t count = 99;
    size_t bytes = 99;
    require_result(leo2_test_balanced_execution_tiles(
            0, 0, &count, &bytes),
        LEO2_SUCCESS, "empty balanced execution geometry");
    require(count == 0 && bytes == 0,
        "empty balanced execution geometry is not empty");
    require_result(leo2_test_balanced_execution_tiles(
            64, 0, &count, &bytes),
        LEO2_INVALID_ARGUMENT, "zero requested execution tile");
    require_result(leo2_test_balanced_execution_tiles(
            65, 64, &count, &bytes),
        LEO2_INVALID_ARGUMENT, "unaligned execution byte count");
    require_result(leo2_test_balanced_execution_tiles(
            128, 65, &count, &bytes),
        LEO2_INVALID_ARGUMENT, "unaligned requested execution tile");
    require_result(leo2_test_balanced_execution_tiles(
            128, 64, NULL, &bytes),
        LEO2_INVALID_ARGUMENT, "null execution tile-count output");
}

class Context
{
public:
    explicit Context(leo2_backend requested)
        : value_(NULL), result_(LEO2_INTERNAL_ERROR)
    {
        leo2_context_options options;
        std::memset(&options, 0, sizeof(options));
        options.struct_size = sizeof(options);
        options.backend = requested;
        options.thread_count = 1;
        result_ = leo2_context_create(&options, &value_);
    }

    ~Context() { leo2_context_destroy(value_); }
    leo2_context* get() const { return value_; }
    leo2_result result() const { return result_; }

private:
    Context(const Context&);
    Context& operator=(const Context&);
    leo2_context* value_;
    leo2_result result_;
};

class Codec
{
public:
    Codec(
        leo2_context* context,
        uint32_t k,
        uint32_t r,
        leo2_profile profile = LEO2_PROFILE_LEGACY_HIGH_V1,
        leo2_field field = LEO2_FIELD_GF16,
        uint32_t flags = 0)
        : value_(NULL)
    {
        leo2_codec_options options;
        std::memset(&options, 0, sizeof(options));
        options.struct_size = sizeof(options);
        options.flags = flags;
        require_result(leo2_codec_create(context, k, r,
            profile, field, flags == 0 ? NULL : &options, &value_),
            LEO2_SUCCESS, "codec create");
    }

    ~Codec() { leo2_codec_destroy(value_); }
    leo2_codec* get() const { return value_; }

private:
    Codec(const Codec&);
    Codec& operator=(const Codec&);
    leo2_codec* value_;
};

class Plan
{
public:
    Plan(leo2_codec* codec, uint32_t k, uint32_t r, uint32_t losses)
        : value_(NULL)
    {
        std::vector<uint8_t> original_present(k, 1);
        std::vector<uint8_t> recovery_present(r, 1);
        for (uint32_t i = 0; i < losses; ++i)
            original_present[i] = 0;
        require_result(leo2_decode_plan_create(
                codec, original_present.data(), recovery_present.data(),
                &value_),
            LEO2_SUCCESS, "decode plan create");
    }

    ~Plan() { leo2_decode_plan_destroy(value_); }
    leo2_decode_plan* get() const { return value_; }

private:
    Plan(const Plan&);
    Plan& operator=(const Plan&);
    leo2_decode_plan* value_;
};

struct ExecutionTiles
{
    size_t count;
    size_t maximum_bytes;
};

ExecutionTiles encode_tiles(const leo2_codec* codec, uint64_t shard_bytes)
{
    ExecutionTiles tiles = { 0, 0 };
    require(leopard::backend::TestOnlyGetEncodeExecutionTiles(
            codec, shard_bytes, &tiles.count, &tiles.maximum_bytes),
        "encode execution-tile query");
    return tiles;
}

ExecutionTiles decode_tiles(
    leo2_codec* codec,
    uint32_t k,
    uint32_t r,
    uint32_t losses,
    uint64_t shard_bytes)
{
    Plan plan(codec, k, r, losses);
    ExecutionTiles tiles = { 0, 0 };
    require(leopard::backend::TestOnlyGetDecodeExecutionTiles(
            plan.get(), shard_bytes, &tiles.count, &tiles.maximum_bytes),
        "decode execution-tile query");
    return tiles;
}

void require_tiles(
    const ExecutionTiles& actual,
    size_t expected_count,
    size_t expected_maximum_bytes,
    const char* message)
{
    if (actual.count != expected_count ||
        actual.maximum_bytes != expected_maximum_bytes)
    {
        throw std::runtime_error(std::string(message) +
            ": count=" + std::to_string(actual.count) +
            " bytes=" + std::to_string(actual.maximum_bytes));
    }
}

void require_cache_policy(
    leo2_context* context,
    uint64_t detected_l3_bytes,
    uint64_t expected_l3_bytes,
    uint64_t expected_target_bytes,
    uint64_t expected_threshold_bytes)
{
    leopard::backend::TestOnlySetContextL3Bytes(
        context, detected_l3_bytes);
    uint64_t l3_bytes = 0;
    uint64_t target_bytes = 0;
    uint64_t threshold_bytes = 0;
    require(leopard::backend::TestOnlyGetContextGF16CachePolicy(
            context, &l3_bytes, &target_bytes, &threshold_bytes),
        "GF16 cache-policy query");
    require(l3_bytes == expected_l3_bytes &&
            target_bytes == expected_target_bytes &&
            threshold_bytes == expected_threshold_bytes,
        "GF16 cache-policy derivation");
}

void require_context_cache_policy_wiring(
    Context& context,
    uint64_t detected_l3_bytes,
    bool cache_sensitive)
{
    if (context.result() != LEO2_SUCCESS)
        return;

    static const uint64_t mib = UINT64_C(1024) * 1024;
    const bool has_gf16 =
        (leo2_context_field_mask(context.get()) &
            LEO2_FIELD_MASK_GF16) != 0;
    const uint64_t detected =
        cache_sensitive && has_gf16 ? detected_l3_bytes : 0;
    uint64_t expected_l3 = 32 * mib;
    if (detected >= expected_l3)
    {
        expected_l3 = detected < 96 * mib ? detected : 96 * mib;
        expected_l3 = expected_l3 / mib * mib;
    }

    uint64_t actual_l3 = 0;
    uint64_t actual_target = 0;
    uint64_t actual_threshold = 0;
    require(leopard::backend::TestOnlyGetContextGF16CachePolicy(
            context.get(), &actual_l3, &actual_target, &actual_threshold),
        "context-created GF16 cache-policy query");
    require(actual_l3 == expected_l3 &&
            actual_target == 16 * mib &&
            actual_threshold ==
                (expected_l3 < 64 * mib ? 64 * mib : expected_l3),
        "context creation did not wire the detected GF16 cache policy");
}

void test_gf16_cache_policy(Context& avx2)
{
    if (avx2.result() != LEO2_SUCCESS ||
        (leo2_context_field_mask(avx2.get()) &
            LEO2_FIELD_MASK_GF16) == 0)
        return;

    static const uint64_t mib = UINT64_C(1024) * 1024;
    uint64_t original_l3_bytes = 0;
    uint64_t original_target_bytes = 0;
    uint64_t original_threshold_bytes = 0;
    require(leopard::backend::TestOnlyGetContextGF16CachePolicy(
            avx2.get(), &original_l3_bytes, &original_target_bytes,
            &original_threshold_bytes),
        "original GF16 cache-policy query");
    (void)original_target_bytes;
    (void)original_threshold_bytes;
    require_cache_policy(avx2.get(), 0, 32 * mib, 16 * mib, 64 * mib);
    require_cache_policy(avx2.get(), 8 * mib, 32 * mib, 16 * mib, 64 * mib);
    require_cache_policy(avx2.get(), 32 * mib, 32 * mib, 16 * mib, 64 * mib);
    require_cache_policy(avx2.get(), 64 * mib, 64 * mib, 16 * mib, 64 * mib);
    require_cache_policy(avx2.get(), 96 * mib, 96 * mib, 16 * mib, 96 * mib);
    require_cache_policy(avx2.get(), 256 * mib, 96 * mib, 16 * mib, 96 * mib);

    const uint32_t forced_specialized =
        LEO2_CODEC_FORCE_SPECIALIZED_DECODE;

    // The 32-MiB fallback exactly reproduces the established generic boundary.
    leopard::backend::TestOnlySetContextL3Bytes(avx2.get(), 32 * mib);
    {
        Codec low(avx2.get(), 64, 193, LEO2_PROFILE_LOW_V1,
            LEO2_FIELD_GF16, forced_specialized);
        require_tiles(encode_tiles(low.get(), 512U * 1024U - 64U),
            1, 512U * 1024U - 64U,
            "32-MiB low encode below threshold");
        require_tiles(encode_tiles(low.get(), 512U * 1024U),
            4, 128U * 1024U,
            "32-MiB low encode at threshold");
        require_tiles(decode_tiles(
                low.get(), 64, 193, 9, 512U * 1024U),
            4, 128U * 1024U,
            "32-MiB low decode at threshold");
    }
    {
        Codec side256(avx2.get(), 1000, 200,
            LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF16,
            forced_specialized);
        require_tiles(encode_tiles(side256.get(), 64U * 1024U),
            2, 32U * 1024U,
            "32-MiB high side-256 encode override");
        require_tiles(decode_tiles(
                side256.get(), 1000, 200, 9, 64U * 1024U),
            4, 16U * 1024U,
            "32-MiB high side-256 decode override");
    }
    {
        Codec side512(avx2.get(), 2000, 500,
            LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF16,
            forced_specialized);
        require_tiles(encode_tiles(side512.get(), 64U * 1024U),
            8, 8U * 1024U,
            "32-MiB high side-512 encode override");
        require_tiles(decode_tiles(
                side512.get(), 2000, 500, 9, 64U * 1024U),
            4, 16U * 1024U,
            "32-MiB high side-512 decode override");
    }
    {
        Codec low_side256(avx2.get(), 200, 800,
            LEO2_PROFILE_LOW_V1, LEO2_FIELD_GF16,
            forced_specialized);
        require_tiles(encode_tiles(low_side256.get(), 64U * 1024U),
            1, 64U * 1024U,
            "32-MiB low side-256 encode inherited a high-only override");
        require_tiles(decode_tiles(
                low_side256.get(), 200, 800, 9, 64U * 1024U),
            1, 64U * 1024U,
            "32-MiB low side-256 unqualified pattern");
        require_tiles(decode_tiles(
                low_side256.get(), 200, 800, 8, 64U * 1024U),
            1, 64U * 1024U,
            "32-MiB low side-256 loss neighbor exclusion");
        require_tiles(decode_tiles(
                low_side256.get(), 200, 800, 9, 96U * 1024U),
            1, 96U * 1024U,
            "32-MiB low side-256 byte neighbor exclusion");
        require_tiles(decode_tiles(
                low_side256.get(), 200, 800, 9, 64U * 1024U + 2U),
            1, 64U * 1024U,
            "32-MiB low side-256 ragged neighbor exclusion");
        require_tiles(decode_tiles(
                low_side256.get(), 200, 800, 9, 128U * 1024U),
            4, 32U * 1024U,
            "32-MiB low side-256 generic crossover");
    }
    {
        Codec low_side256_neighbor(avx2.get(), 201, 799,
            LEO2_PROFILE_LOW_V1, LEO2_FIELD_GF16,
            forced_specialized);
        require_tiles(decode_tiles(
                low_side256_neighbor.get(), 201, 799, 9, 64U * 1024U),
            1, 64U * 1024U,
            "32-MiB low side-256 count neighbor exclusion");
    }

    // A 96-MiB context declines the measured 32/64-MiB over-tiling region.
    leopard::backend::TestOnlySetContextL3Bytes(avx2.get(), 96 * mib);
    {
        Codec side256(avx2.get(), 1000, 200,
            LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF16,
            forced_specialized);
        require_tiles(encode_tiles(side256.get(), 64U * 1024U),
            1, 64U * 1024U,
            "96-MiB high side-256 encode exclusion");
        require_tiles(decode_tiles(
                side256.get(), 1000, 200, 9, 64U * 1024U),
            1, 64U * 1024U,
            "96-MiB high side-256 decode exclusion");
    }
    {
        Codec side512(avx2.get(), 2000, 500,
            LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF16,
            forced_specialized);
        require_tiles(encode_tiles(side512.get(), 64U * 1024U),
            1, 64U * 1024U,
            "96-MiB high side-512 encode exclusion");
        require_tiles(decode_tiles(
                side512.get(), 2000, 500, 9, 64U * 1024U),
            1, 64U * 1024U,
            "96-MiB high side-512 decode exclusion");
        require_tiles(encode_tiles(side512.get(), 128U * 1024U),
            6, 21888,
            "96-MiB high side-512 scaled encode override");
        require_tiles(decode_tiles(
                side512.get(), 2000, 500, 9, 128U * 1024U),
            3, 43712,
            "96-MiB high side-512 scaled decode override");
    }
    {
        Codec side512_max_loss(avx2.get(), 2000, 512,
            LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF16,
            forced_specialized);
        require_tiles(decode_tiles(
                side512_max_loss.get(), 2000, 512, 512, 64U * 1024U),
            2, 32U * 1024U,
            "96-MiB high side-512 maximum-loss live-set inclusion");
    }
    {
        /*
            A codec-level one-shot query reserves 2T+R rows, while this
            small-loss plan retains only 2T+9.  At 376 KiB the codec query's
            minimum 2T policy rows remain below the 96-MiB threshold while
            the known plan crosses it.  The codec query must therefore keep a
            conservative full pass; otherwise leo2_decode() could build the
            tiled plan and reject its caller-provided scratch.
        */
        Codec one_shot_bound(avx2.get(), 300, 100,
            LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF16,
            forced_specialized);
        Plan small_loss(
            one_shot_bound.get(), 300, 100, 9);
        static const size_t bytes = 376U * 1024U;
        size_t codec_scratch_bytes = 0;
        size_t plan_scratch_bytes = 0;
        require_result(leo2_decode_scratch_size(
                one_shot_bound.get(), bytes, &codec_scratch_bytes),
            LEO2_SUCCESS, "96-MiB codec decode scratch query");
        require_result(leo2_decode_plan_scratch_size(
                small_loss.get(), bytes, &plan_scratch_bytes),
            LEO2_SUCCESS, "96-MiB plan decode scratch query");
        require_tiles(decode_tiles(
                one_shot_bound.get(), 300, 100, 9, bytes),
            7, 55040, "96-MiB small-loss tiled decode");
        require(codec_scratch_bytes >= plan_scratch_bytes,
            "codec decode scratch query under-sized a small-loss plan");
    }
    {
        Codec low(avx2.get(), 64, 193, LEO2_PROFILE_LOW_V1,
            LEO2_FIELD_GF16, forced_specialized);
        require_tiles(encode_tiles(low.get(), 768U * 1024U - 64U),
            1, 768U * 1024U - 64U,
            "96-MiB low encode below threshold");
        require_tiles(encode_tiles(low.get(), 768U * 1024U),
            6, 128U * 1024U,
            "96-MiB low encode at threshold");
        require_tiles(decode_tiles(
                low.get(), 64, 193, 9, 768U * 1024U),
            6, 128U * 1024U,
            "96-MiB low decode at threshold");
    }
    {
        Codec low_side512(avx2.get(), 300, 800,
            LEO2_PROFILE_LOW_V1, LEO2_FIELD_GF16,
            forced_specialized);
        require_tiles(encode_tiles(low_side512.get(), 96U * 1024U),
            6, 16U * 1024U,
            "96-MiB low side-512 generic encode crossover");
        require_tiles(decode_tiles(
                low_side512.get(), 300, 800, 9, 96U * 1024U),
            6, 16U * 1024U,
            "96-MiB low side-512 generic decode crossover");
        require_tiles(encode_tiles(
                low_side512.get(), 128U * 1024U),
            8, 16U * 1024U,
            "96-MiB low side-512 generic encode threshold");
        require_tiles(decode_tiles(
                low_side512.get(), 300, 800, 9, 128U * 1024U),
            8, 16U * 1024U,
            "96-MiB low side-512 generic decode threshold");
    }
    {
        Codec low_side256(avx2.get(), 200, 800,
            LEO2_PROFILE_LOW_V1, LEO2_FIELD_GF16,
            forced_specialized);
        require_tiles(decode_tiles(
                low_side256.get(), 200, 800, 9, 64U * 1024U),
            1, 64U * 1024U,
            "96-MiB low side-256 cache exclusion");
    }

#ifdef LEO_HAS_FF8
    // GF8 owns a separately calibrated table and must ignore the GF16 budget.
    leopard::backend::TestOnlySetContextL3Bytes(avx2.get(), 32 * mib);
    Codec gf8_32(avx2.get(), 100, 128, LEO2_PROFILE_LOW_V1,
        LEO2_FIELD_GF8, forced_specialized);
    const ExecutionTiles gf8_encode_32 =
        encode_tiles(gf8_32.get(), 512U * 1024U);
    const ExecutionTiles gf8_decode_32 = decode_tiles(
        gf8_32.get(), 100, 128, 9, 512U * 1024U);
    leopard::backend::TestOnlySetContextL3Bytes(avx2.get(), 96 * mib);
    require_tiles(encode_tiles(gf8_32.get(), 512U * 1024U),
        gf8_encode_32.count, gf8_encode_32.maximum_bytes,
        "GF8 encode changed with GF16 cache budget");
    require_tiles(decode_tiles(
            gf8_32.get(), 100, 128, 9, 512U * 1024U),
        gf8_decode_32.count, gf8_decode_32.maximum_bytes,
        "GF8 decode changed with GF16 cache budget");
#endif
    leopard::backend::TestOnlySetContextL3Bytes(
        avx2.get(), original_l3_bytes);
}

#if defined(__linux__)
class TemporaryCacheTree
{
public:
    TemporaryCacheTree()
    {
        char pattern[] = "/tmp/leopard2-cache-topology-XXXXXX";
        char* created = mkdtemp(pattern);
        require(created != NULL, "temporary cache root");
        root_ = created;
        directories_.push_back(root_);
    }

    ~TemporaryCacheTree()
    {
        for (std::vector<std::string>::reverse_iterator i = files_.rbegin();
             i != files_.rend(); ++i)
            unlink(i->c_str());
        for (std::vector<std::string>::reverse_iterator i =
                 directories_.rbegin(); i != directories_.rend(); ++i)
            rmdir(i->c_str());
    }

    const std::string& root() const { return root_; }

    void directory(const std::string& relative)
    {
        const std::string path = root_ + "/" + relative;
        require(mkdir(path.c_str(), 0700) == 0,
            "temporary cache directory");
        directories_.push_back(path);
    }

    void write(const std::string& relative, const char* text)
    {
        const std::string path = root_ + "/" + relative;
        FILE* file = std::fopen(path.c_str(), "wb");
        require(file != NULL, "temporary cache attribute create");
        const size_t bytes = std::strlen(text);
        const bool wrote = std::fwrite(text, 1, bytes, file) == bytes;
        const bool closed = std::fclose(file) == 0;
        require(wrote && closed, "temporary cache attribute write");
        files_.push_back(path);
    }

private:
    std::string root_;
    std::vector<std::string> directories_;
    std::vector<std::string> files_;
};

void add_fake_cache_cpu(
    TemporaryCacheTree& tree,
    uint32_t cpu,
    const char* shared_l3_cpus,
    const char* l3_size)
{
    const std::string prefix = "cpu" + std::to_string(cpu);
    tree.directory(prefix);
    tree.directory(prefix + "/cache");
    for (unsigned index = 0; index < 4; ++index)
        tree.directory(prefix + "/cache/index" + std::to_string(index));

    tree.write(prefix + "/cache/index0/type", "Data\n");
    tree.write(prefix + "/cache/index0/level", "1\n");
    tree.write(prefix + "/cache/index0/size", "32K\n");
    tree.write(prefix + "/cache/index0/shared_cpu_list",
        (std::to_string(cpu) + "\n").c_str());

    tree.write(prefix + "/cache/index1/type", "Instruction\n");

    tree.write(prefix + "/cache/index2/type", "Unified\n");
    tree.write(prefix + "/cache/index2/level", "2\n");
    tree.write(prefix + "/cache/index2/size", "1024K\n");
    tree.write(prefix + "/cache/index2/shared_cpu_list",
        (std::to_string(cpu) + "\n").c_str());

    tree.write(prefix + "/cache/index3/type", "Unified\n");
    tree.write(prefix + "/cache/index3/level", "3\n");
    tree.write(prefix + "/cache/index3/size", l3_size);
    tree.write(prefix + "/cache/index3/shared_cpu_list", shared_l3_cpus);
}

void test_linux_cache_topology()
{
    struct ParseCase
    {
        const char* text;
        bool valid;
        uint64_t bytes;
    };
    static const ParseCase cases[] = {
        { "32K\n", true, UINT64_C(32) * 1024 },
        { "98304K\n", true, UINT64_C(96) * 1024 * 1024 },
        { "96M", true, UINT64_C(96) * 1024 * 1024 },
        { "1G\n", true, UINT64_C(1) << 30 },
        { "1048576\n", true, UINT64_C(1) << 20 },
        { "", false, 0 },
        { "0K\n", false, 0 },
        { "-1M\n", false, 0 },
        { "+1M\n", false, 0 },
        { "1.5M\n", false, 0 },
        { "96m\n", false, 0 },
        { " 96M\n", false, 0 },
        { "96M \n", false, 0 },
        { "96M\n\n", false, 0 },
        { "96M\r\n", false, 0 },
        { "18446744073709551615G\n", false, 0 }
    };
    for (size_t i = 0; i < sizeof(cases) / sizeof(cases[0]); ++i)
    {
        uint64_t bytes = UINT64_MAX;
        const bool valid = leopard::backend::TestOnlyParseLinuxCacheSize(
            cases[i].text, &bytes);
        require(valid == cases[i].valid,
            "Linux cache-size parser validity");
        if (valid)
            require(bytes == cases[i].bytes,
                "Linux cache-size parser value");
    }

    TemporaryCacheTree tree;
    add_fake_cache_cpu(tree, 0, "0-1\n", "98304K\n");
    add_fake_cache_cpu(tree, 1, "0-1\n", "96M\n");
    add_fake_cache_cpu(tree, 8, "8\n", "32768K\n");
    tree.directory("cpu8/cache/index4");
    tree.write("cpu8/cache/index4/type", "Unified\n");
    tree.write("cpu8/cache/index4/level", "4\n");
    tree.write("cpu8/cache/index4/size", "256M\n");
    tree.write("cpu8/cache/index4/shared_cpu_list", "8\n");
    add_fake_cache_cpu(tree, 9, "9\n", "1K\n");

    const uint32_t large_domain[] = { 0, 1 };
    uint64_t bytes = 0;
    require(leopard::backend::TestOnlyDetectLinuxL3Bytes(
            tree.root().c_str(), large_domain, 2, &bytes) &&
            bytes == UINT64_C(96) * 1024 * 1024,
        "large Linux cache-domain detection");

    const uint32_t mixed_domains[] = { 0, 1, 8 };
    require(leopard::backend::TestOnlyDetectLinuxL3Bytes(
            tree.root().c_str(), mixed_domains, 3, &bytes) &&
            bytes == UINT64_C(32) * 1024 * 1024,
        "mixed Linux cache-domain minimum");
    const uint32_t l4_domain[] = { 8 };
    require(leopard::backend::TestOnlyDetectLinuxL3Bytes(
            tree.root().c_str(), l4_domain, 1, &bytes) &&
            bytes == UINT64_C(32) * 1024 * 1024,
        "Linux L4 replaced calibrated L3 capacity");
    const uint32_t implausible_l3[] = { 9 };
    require(!leopard::backend::TestOnlyDetectLinuxL3Bytes(
            tree.root().c_str(), implausible_l3, 1, &bytes),
        "implausibly small Linux L3 accepted");

    const uint32_t unsorted[] = { 8, 1 };
    require(!leopard::backend::TestOnlyDetectLinuxL3Bytes(
            tree.root().c_str(), unsorted, 2, &bytes),
        "unsorted Linux CPU list accepted");
    require(!leopard::backend::TestOnlyDetectLinuxL3Bytes(
            tree.root().c_str(), mixed_domains, 3, NULL),
        "null Linux cache output accepted");
}
#else
void test_linux_cache_topology()
{
    uint64_t bytes = UINT64_C(0x123456789abcdef0);
    require(!leopard::backend::TestOnlyParseLinuxCacheSize(
            "32M\n", &bytes) &&
            bytes == UINT64_C(0x123456789abcdef0),
        "non-Linux cache-size hook modified its output");
    const uint32_t cpu = 0;
    require(!leopard::backend::TestOnlyDetectLinuxL3Bytes(
            "/non-linux", &cpu, 1, &bytes) &&
            bytes == UINT64_C(0x123456789abcdef0),
        "non-Linux cache-topology hook modified its output");
}
#endif

class AlignedBuffer
{
public:
    explicit AlignedBuffer(size_t bytes)
        : value_(NULL)
    {
#if defined(_MSC_VER)
        value_ = _aligned_malloc(bytes, leo2_scratch_alignment());
#else
        if (posix_memalign(&value_, leo2_scratch_alignment(), bytes) != 0)
            value_ = NULL;
#endif
        if (!value_)
            throw std::bad_alloc();
        std::memset(value_, 0, bytes);
    }

    ~AlignedBuffer()
    {
#if defined(_MSC_VER)
        _aligned_free(value_);
#else
        std::free(value_);
#endif
    }

    void* get() const { return value_; }

private:
    AlignedBuffer(const AlignedBuffer&);
    AlignedBuffer& operator=(const AlignedBuffer&);
    void* value_;
};

leo2_backend selected_backend(
    const leo2_codec* codec,
    uint64_t bytes,
    uint32_t count,
    uint32_t prefix)
{
    leo2_backend backend = LEO2_BACKEND_AUTO;
    require_result(leo2_test_codec_transform_encode_backend(
        codec, bytes, count, prefix, &backend), LEO2_SUCCESS,
        "encode backend query");
    return backend;
}

Shards make_original(uint32_t k, size_t bytes)
{
    Shards original(k, Bytes(bytes));
    uint32_t state = 0x243f6a88U;
    for (uint32_t shard = 0; shard < k; ++shard)
        for (size_t i = 0; i < bytes; ++i)
        {
            state = state * 1664525U + 1013904223U;
            original[shard][i] = static_cast<uint8_t>(
                (state >> 24) ^ shard ^ static_cast<uint32_t>(i));
        }
    return original;
}

Shards encode(
    const leo2_codec* codec,
    const Shards& original,
    uint32_t r,
    bool sparse)
{
    const size_t bytes = original[0].size();
    std::vector<const void*> original_ptrs(original.size());
    for (size_t i = 0; i < original.size(); ++i)
        original_ptrs[i] = original[i].data();
    Shards recovery(r, Bytes(bytes));
    std::vector<void*> recovery_ptrs(r);
    for (uint32_t i = 0; i < r; ++i)
        recovery_ptrs[i] = sparse && (i & 1U) != 0
            ? NULL : recovery[i].data();
    size_t scratch_bytes = 0;
    require_result(leo2_encode_scratch_size(codec, bytes, &scratch_bytes),
        LEO2_SUCCESS, "scratch query");
    AlignedBuffer scratch(scratch_bytes);
    require_result(leo2_encode(codec, bytes, original_ptrs.data(),
        recovery_ptrs.data(), scratch.get(), scratch_bytes),
        LEO2_SUCCESS, "encode");
    return recovery;
}

void run_t8_mixed_binding(
    leo2_codec* codec,
    uint32_t k,
    uint32_t r,
    const size_t* byte_counts,
    size_t item_count,
    uint64_t expected_whole_transform_calls,
    uint64_t expected_ifft_groups_per_call)
{
    std::vector<Shards> original(item_count);
    std::vector<Shards> expected(item_count);
    std::vector<Shards> actual(item_count);
    std::vector<std::vector<const void*> > original_ptrs(
        item_count, std::vector<const void*>(k, NULL));
    std::vector<std::vector<void*> > recovery_ptrs(
        item_count, std::vector<void*>(r, NULL));
    std::vector<std::unique_ptr<AlignedBuffer> > scratch(item_count);
    std::vector<leo2_encode_batch_item> items(item_count);

    for (size_t item_index = 0; item_index < item_count; ++item_index)
    {
        const size_t bytes = byte_counts[item_index];
        original[item_index] = make_original(k, bytes);
        expected[item_index] = encode(
            codec, original[item_index], r, false);
        actual[item_index] = Shards(r, Bytes(bytes, 0xa5));
        for (uint32_t source = 0; source < k; ++source)
            original_ptrs[item_index][source] =
                original[item_index][source].data();
        for (uint32_t parity = 0; parity < r; ++parity)
            recovery_ptrs[item_index][parity] =
                actual[item_index][parity].data();
        size_t scratch_bytes = 0;
        require_result(leo2_encode_scratch_size(
                codec, bytes, &scratch_bytes),
            LEO2_SUCCESS, "mixed T8 scratch query");
        scratch[item_index].reset(new AlignedBuffer(scratch_bytes));
        items[item_index].shard_bytes = bytes;
        items[item_index].original = original_ptrs[item_index].data();
        items[item_index].recovery = recovery_ptrs[item_index].data();
        items[item_index].scratch = scratch[item_index]->get();
        items[item_index].scratch_bytes = scratch_bytes;
    }

    leo2_encode_batch_binding* binding = NULL;
    require_result(leo2_encode_batch_binding_create(
            codec, items.data(), items.size(), &binding),
        LEO2_SUCCESS, "mixed T8 binding create");
    leopard::ff8::TestOnlyResetHighEncodeCounts();
    require_result(leo2_encode_batch_binding_execute(binding),
        LEO2_SUCCESS, "mixed T8 binding execute");
    const leopard::ff8::TestOnlyHighEncodeCounts counts =
        leopard::ff8::TestOnlyGetHighEncodeCounts();
    if (counts.whole_transform_calls != expected_whole_transform_calls ||
        counts.forward_fused_calls != expected_whole_transform_calls ||
        (expected_whole_transform_calls != 0 &&
         counts.ifft_butterfly4_out_of_place !=
            expected_ifft_groups_per_call *
                expected_whole_transform_calls))
    {
        throw std::runtime_error(
            "mixed T8 binding executed the wrong transform route: whole=" +
            std::to_string(counts.whole_transform_calls) +
            " forward=" + std::to_string(counts.forward_fused_calls) +
            " ifft4=" +
            std::to_string(counts.ifft_butterfly4_out_of_place) +
            " expected_whole=" +
            std::to_string(expected_whole_transform_calls));
    }
    require(actual == expected, "mixed T8 binding changed parity bytes");
    leo2_encode_batch_binding_destroy(binding);
}

void test_t8_one_block_mixed_binding(Context& avx2)
{
    if (avx2.result() != LEO2_SUCCESS)
        return;
    Codec codec(avx2.get(), 5, 5,
        LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8);
    static const size_t qualified[] = { 64, 128, 320, 512 };
    leopard2_internal::CodecEncodePathInfo path = {};
    require(leopard2_internal::GetCodecEncodePathInfo(
            codec.get(), 512, 5, &path),
        "mixed T8 path query");
    const uint64_t expected_qualified_calls =
        path.high_t8_partial_binding_selected ? 4 : 0;
    run_t8_mixed_binding(
        codec.get(), 5, 5, qualified,
        sizeof(qualified) / sizeof(qualified[0]),
        expected_qualified_calls, 2);

    /*
        Exercise the fixed-shape K=5,R=5 kernel in one heterogeneous binding
        with the original 64-byte tile.  Both items must select the same
        immutable AVX2 callback while dispatching their byte-specific bodies.
    */
    static const size_t exact_kernel[] = { 64, 1088 };
    path = leopard2_internal::CodecEncodePathInfo();
    require(leopard2_internal::GetCodecEncodePathInfo(
            codec.get(), 1088, 5, &path),
        "mixed exact T8 path query");
    const uint64_t expected_exact_calls =
        path.high_t8_partial_binding_selected ? 2 : 0;
    run_t8_mixed_binding(
        codec.get(), 5, 5, exact_kernel,
        sizeof(exact_kernel) / sizeof(exact_kernel[0]),
        expected_exact_calls, 2);

    /*
        One unqualified item must disable the binding-level shortcut for the
        complete heterogeneous batch.  The mature per-item path still emits
        identical parity for both items.
    */
    static const size_t boundary[] = { 64, 513 };
    path = leopard2_internal::CodecEncodePathInfo();
    require(leopard2_internal::GetCodecEncodePathInfo(
            codec.get(), 513, 5, &path) &&
            !path.high_t8_partial_binding_selected,
        "mixed T8 boundary unexpectedly selected the shortcut");
    run_t8_mixed_binding(
        codec.get(), 5, 5, boundary,
        sizeof(boundary) / sizeof(boundary[0]), 0, 2);
}

void test_t8_two_block_mixed_binding(Context& avx2)
{
    if (avx2.result() != LEO2_SUCCESS)
        return;
    Codec codec(avx2.get(), 13, 5,
        LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8);
    static const size_t qualified[] = { 64, 128, 320, 512 };
    leopard2_internal::CodecEncodePathInfo path = {};
    require(leopard2_internal::GetCodecEncodePathInfo(
            codec.get(), 512, 5, &path),
        "mixed two-block T8 path query");
    const uint64_t expected_qualified_calls =
        path.high_t8_two_block_binding_selected ? 4 : 0;
    run_t8_mixed_binding(
        codec.get(), 13, 5, qualified,
        sizeof(qualified) / sizeof(qualified[0]),
        expected_qualified_calls, 4);

    /*
        As with the one-block route, a single unqualified item must disable
        the binding-level callback for the whole heterogeneous batch.
    */
    static const size_t boundary[] = { 64, 513 };
    run_t8_mixed_binding(
        codec.get(), 13, 5, boundary,
        sizeof(boundary) / sizeof(boundary[0]), 0, 4);
}

Shards encode_unaligned(
    const leo2_codec* codec,
    const Shards& original,
    uint32_t r)
{
    const size_t bytes = original[0].size();
    Shards original_storage(original.size(), Bytes(bytes + 3));
    std::vector<const void*> original_ptrs(original.size());
    for (size_t i = 0; i < original.size(); ++i)
    {
        std::memcpy(original_storage[i].data() + 1,
            original[i].data(), bytes);
        original_ptrs[i] = original_storage[i].data() + 1;
    }

    Shards recovery_storage(r, Bytes(bytes + 5));
    std::vector<void*> recovery_ptrs(r);
    for (uint32_t i = 0; i < r; ++i)
        recovery_ptrs[i] = recovery_storage[i].data() + 3;

    size_t scratch_bytes = 0;
    require_result(leo2_encode_scratch_size(codec, bytes, &scratch_bytes),
        LEO2_SUCCESS, "unaligned scratch query");
    AlignedBuffer scratch(scratch_bytes);
    require_result(leo2_encode(codec, bytes, original_ptrs.data(),
        recovery_ptrs.data(), scratch.get(), scratch_bytes),
        LEO2_SUCCESS, "unaligned encode");

    Shards recovery(r, Bytes(bytes));
    for (uint32_t i = 0; i < r; ++i)
        std::memcpy(recovery[i].data(),
            recovery_storage[i].data() + 3, bytes);
    return recovery;
}

void require_encode_overlap_rejected(
    const leo2_codec* codec,
    Shards& original,
    uint32_t r)
{
    const size_t bytes = original[0].size();
    const Shards original_before = original;
    std::vector<const void*> original_ptrs(original.size());
    for (size_t i = 0; i < original.size(); ++i)
        original_ptrs[i] = original[i].data();
    Shards recovery(r, Bytes(bytes));
    std::vector<void*> recovery_ptrs(r);
    for (uint32_t i = 0; i < r; ++i)
        recovery_ptrs[i] = recovery[i].data();
    recovery_ptrs[0] = original[0].data();
    size_t scratch_bytes = 0;
    require_result(leo2_encode_scratch_size(codec, bytes, &scratch_bytes),
        LEO2_SUCCESS, "overlap scratch query");
    AlignedBuffer scratch(scratch_bytes);
    require_result(leo2_encode(codec, bytes, original_ptrs.data(),
        recovery_ptrs.data(), scratch.get(), scratch_bytes),
        LEO2_OVERLAP, "encode input/output overlap");
    require(original == original_before,
        "rejected encode overlap modified a source shard");
}

Shards encode_prefix(
    const leo2_codec* codec,
    const Shards& original,
    uint32_t r,
    uint32_t prefix)
{
    const size_t bytes = original[0].size();
    std::vector<const void*> original_ptrs(original.size());
    for (size_t i = 0; i < original.size(); ++i)
        original_ptrs[i] = original[i].data();
    Shards recovery(r, Bytes(bytes));
    std::vector<void*> recovery_ptrs(r, NULL);
    for (uint32_t i = 0; i < prefix; ++i)
        recovery_ptrs[i] = recovery[i].data();
    size_t scratch_bytes = 0;
    require_result(leo2_encode_scratch_size(codec, bytes, &scratch_bytes),
        LEO2_SUCCESS, "prefix scratch query");
    AlignedBuffer scratch(scratch_bytes);
    require_result(leo2_encode(codec, bytes, original_ptrs.data(),
        recovery_ptrs.data(), scratch.get(), scratch_bytes),
        LEO2_SUCCESS, "prefix encode");
    return recovery;
}

Shards encode_legacy(const Shards& original, uint32_t r)
{
    const size_t bytes = original[0].size();
    const uint32_t k = static_cast<uint32_t>(original.size());
    std::vector<const void*> original_ptrs(k);
    for (uint32_t i = 0; i < k; ++i)
        original_ptrs[i] = original[i].data();
    const unsigned work_count = leo_encode_work_count(k, r);
    require(work_count >= r, "legacy encode work count");
    Shards work(work_count, Bytes(bytes));
    std::vector<void*> work_ptrs(work_count);
    for (unsigned i = 0; i < work_count; ++i)
        work_ptrs[i] = work[i].data();
    require(leo_encode(bytes, k, r, work_count,
                original_ptrs.data(), work_ptrs.data()) == Leopard_Success,
        "legacy encode");
    work.resize(r);
    return work;
}

void require_sparse_matches_full(
    const Shards& sparse, const Shards& full)
{
    require(sparse.size() == full.size(), "sparse parity size mismatch");
    for (size_t i = 0; i < full.size(); i += 2)
        require(sparse[i] == full[i], "sparse parity differs from full parity");
}

void run_injected_gf16_cache_round_trip(
    Context& context,
    uint32_t k,
    uint32_t r,
    leo2_profile profile,
    size_t shard_bytes)
{
    static const uint32_t losses = 9;
    Codec codec(context.get(), k, r, profile, LEO2_FIELD_GF16,
        LEO2_CODEC_FORCE_SPECIALIZED_DECODE);
    const Shards original = make_original(k, shard_bytes);
    static const uint64_t mib = UINT64_C(1024) * 1024;
    leopard::backend::TestOnlySetContextL3Bytes(context.get(), 32 * mib);
    const Shards established_recovery =
        encode(codec.get(), original, r, false);
    leopard::backend::TestOnlySetContextL3Bytes(context.get(), 96 * mib);
    const ExecutionTiles encode_geometry =
        encode_tiles(codec.get(), shard_bytes);
    const ExecutionTiles decode_geometry =
        decode_tiles(codec.get(), k, r, losses, shard_bytes);
    require(encode_geometry.count > 1 &&
            encode_geometry.maximum_bytes < shard_bytes,
        "injected 96-MiB encode did not exercise byte tiling");
    require(decode_geometry.count > 1 &&
            decode_geometry.maximum_bytes < shard_bytes,
        "injected 96-MiB decode did not exercise byte tiling");
    const Shards recovery = encode(codec.get(), original, r, false);
    require(recovery == established_recovery,
        "injected GF16 cache policy changed encoded parity bytes");

    Plan plan(codec.get(), k, r, losses);
    std::vector<const void*> original_inputs(k);
    for (uint32_t i = 0; i < k; ++i)
        original_inputs[i] = i < losses ? NULL : original[i].data();
    std::vector<const void*> recovery_inputs(r);
    for (uint32_t i = 0; i < r; ++i)
        recovery_inputs[i] = recovery[i].data();

    Shards restored(losses, Bytes(shard_bytes));
    std::vector<void*> restored_outputs(k, NULL);
    for (uint32_t i = 0; i < losses; ++i)
        restored_outputs[i] = restored[i].data();

    size_t scratch_bytes = 0;
    require_result(leo2_decode_plan_scratch_size(
            plan.get(), shard_bytes, &scratch_bytes),
        LEO2_SUCCESS, "injected GF16 decode scratch query");
    AlignedBuffer scratch(scratch_bytes);
    require_result(leo2_decode_plan_execute(
            plan.get(), shard_bytes, original_inputs.data(),
            recovery_inputs.data(), restored_outputs.data(),
            scratch.get(), scratch_bytes),
        LEO2_SUCCESS, "injected GF16 decode execute");
    for (uint32_t i = 0; i < losses; ++i)
        require(restored[i] == original[i],
            "injected GF16 cache policy changed recovered bytes");
}

void run_codec_scratch_query_round_trip(Context& context)
{
    static const uint32_t k = 129;
    static const uint32_t r = 100;
    static const uint32_t losses = 9;
    static const size_t shard_bytes = 224U * 1024U;
    static const uint64_t mib = UINT64_C(1024) * 1024;

    leopard::backend::TestOnlySetContextL3Bytes(context.get(), 32 * mib);
    Codec codec(context.get(), k, r, LEO2_PROFILE_LEGACY_HIGH_V1,
        LEO2_FIELD_GF16, LEO2_CODEC_FORCE_SPECIALIZED_DECODE);
    const Shards original = make_original(k, shard_bytes);
    const Shards recovery = encode(codec.get(), original, r, false);

    std::vector<uint8_t> original_present(k, 1);
    std::vector<uint8_t> recovery_present(r, 1);
    std::vector<const void*> original_inputs(k);
    std::vector<const void*> recovery_inputs(r);
    for (uint32_t i = 0; i < k; ++i)
    {
        if (i < losses)
            original_present[i] = 0;
        original_inputs[i] = i < losses ? NULL : original[i].data();
    }
    for (uint32_t i = 0; i < r; ++i)
        recovery_inputs[i] = recovery[i].data();

    Shards restored(losses, Bytes(shard_bytes));
    std::vector<void*> restored_outputs(k, NULL);
    for (uint32_t i = 0; i < losses; ++i)
        restored_outputs[i] = restored[i].data();

    size_t scratch_bytes = 0;
    require_result(leo2_decode_scratch_size(
            codec.get(), shard_bytes, &scratch_bytes),
        LEO2_SUCCESS, "one-shot GF16 decode scratch query");
    AlignedBuffer scratch(scratch_bytes);
    require_result(leo2_decode(
            codec.get(), shard_bytes, original_present.data(),
            recovery_present.data(), original_inputs.data(),
            recovery_inputs.data(), restored_outputs.data(),
            scratch.get(), scratch_bytes),
        LEO2_SUCCESS, "one-shot GF16 decode at queried scratch bound");
    for (uint32_t i = 0; i < losses; ++i)
        require(restored[i] == original[i],
            "one-shot GF16 queried scratch changed recovered bytes");
}

void test_gf16_cache_policy_bytes(Context& avx2)
{
    if (avx2.result() != LEO2_SUCCESS ||
        (leo2_context_field_mask(avx2.get()) &
            LEO2_FIELD_MASK_GF16) == 0)
        return;

    uint64_t original_l3_bytes = 0;
    uint64_t original_target_bytes = 0;
    uint64_t original_threshold_bytes = 0;
    require(leopard::backend::TestOnlyGetContextGF16CachePolicy(
            avx2.get(), &original_l3_bytes, &original_target_bytes,
            &original_threshold_bytes),
        "original GF16 byte-test cache-policy query");
    (void)original_target_bytes;
    (void)original_threshold_bytes;

    /*
        First exercise the public codec-level scratch bound in the 32-MiB
        max-row/min-row crossover.  The high and low profile cases then each
        reach a 128-MiB live set, at 64- and 128-KiB shards respectively.
        Under the injected 96-MiB policy they execute eight balanced passes
        rather than merely inspecting geometry.
    */
    run_codec_scratch_query_round_trip(avx2);
    run_injected_gf16_cache_round_trip(
        avx2, 1025, 513, LEO2_PROFILE_LEGACY_HIGH_V1, 64U * 1024U);
    run_injected_gf16_cache_round_trip(
        avx2, 257, 513, LEO2_PROFILE_LOW_V1, 128U * 1024U);
    leopard::backend::TestOnlySetContextL3Bytes(
        avx2.get(), original_l3_bytes);
}

void require_explicit_backend(Context& context, leo2_backend expected)
{
    if (context.result() != LEO2_SUCCESS)
        return;
    require(leo2_context_backend(context.get()) == expected,
        "explicit context reported the wrong backend");
    Codec codec(context.get(), 1000, 200);
    require(selected_backend(codec.get(), 4U * 1024U * 1024U + 64U,
                200, 200) == expected,
        "explicit backend was changed by the AUTO calibration bounds");
}

void test_small_high_encode(
    Context& scalar,
    Context& ssse3,
    Context& avx2,
    Context& avx512)
{
    if (avx2.result() != LEO2_SUCCESS)
        return;

    struct TestCase
    {
        uint32_t k;
        uint32_t r;
        size_t bytes;
    };
    static const TestCase cases[] = {
        // T=2 crosses at 2 KiB for every valid K.  Cover the minimum K, an
        // odd message count and ragged byte tail, and the largest GF8 K.
        { 2, 2, 2048 },
        { 3, 2, 2049 },
        { 254, 2, 2048 },
        // T=4 uses the evidence-derived K/byte staircase.  The new lower cells
        // fuse only the final forward butterfly; the existing thresholds keep
        // their register-fused or mature coarse implementations.
        { 3, 4, 2048 },
        { 3, 3, 2049 },
        { 3, 4, 4095 },
        { 3, 4, 4096 },
        { 3, 3, 4097 },
        { 4, 4, 2048 },
        { 4, 3, 4097 },
        { 4, 4, 8191 },
        { 4, 4, 8192 },
        { 4, 3, 8193 },
        { 5, 4, 2048 },
        { 6, 3, 2049 },
        { 7, 4, 2048 },
        { 7, 4, 4095 },
        { 7, 4, 4096 },
        { 16, 4, 2048 },
        { 16, 3, 2049 },
        { 8, 4, 2048 },
        { 8, 4, 4095 },
        { 8, 4, 4096 },
        { 8, 3, 4097 },
        { 9, 4, 2048 },
        { 10, 3, 2049 },
        { 11, 4, 2048 },
        // Final T=4 exact-main fallback-map closure cells.
        { 12, 3, 2048 },
        { 12, 4, 2048 },
        { 12, 3, 2049 },
        { 15, 3, 2048 },
        { 15, 4, 2048 },
        // Retain the previously qualified interior and upper-bound cases.
        { 64, 2, 4096 },
        { 64, 3, 4097 },
        { 240, 4, 65536 },
        { 251, 3, 65536 },
        { 252, 4, 65536 }
    };

    for (size_t case_i = 0;
         case_i < sizeof(cases) / sizeof(cases[0]); ++case_i)
    {
        const TestCase& test_case = cases[case_i];
        const Shards original = make_original(test_case.k, test_case.bytes);
        const Shards original_before = original;
        Codec avx2_codec(avx2.get(), test_case.k, test_case.r,
            LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8);

        leopard::ff8::TestOnlyResetHighEncodeCounts();
        const Shards actual = encode(
            avx2_codec.get(), original, test_case.r, false);
        const leopard::ff8::TestOnlyHighEncodeCounts route =
            leopard::ff8::TestOnlyGetHighEncodeCounts();
        // The aligned prefix is one execution pass.  A ragged 64-byte staging
        // pass deliberately remains below the production byte threshold.
        const uint64_t expected_passes = 1;
        const uint32_t side = test_case.r == 2 ? 2U : 4U;
        const bool fused_t4 = side == 4 &&
            (((test_case.k == 5 || test_case.k == 6 ||
                test_case.k >= 9) && test_case.bytes >= 2U * 1024U) ||
             ((test_case.k == 3 || test_case.k == 7) &&
                test_case.bytes >= 4U * 1024U) ||
             (test_case.k == 4 && test_case.bytes >= 8U * 1024U)) &&
            ((test_case.k >= 3 && test_case.k <= 7) ||
             (test_case.k >= 9 && test_case.k <= 11));
        const uint64_t expected_input_copies =
            (fused_t4 || test_case.k <= side ||
                    (side == 4 && test_case.k % side == 3)
                ? 0U : test_case.k % side) * expected_passes +
            ((test_case.bytes & 63U) != 0 ? test_case.k : 0U);
        require(route.small_transform_calls == expected_passes,
            "dense T=2/T=4 AVX2 encode missed the coarse kernel");
        require(route.input_copy_shards ==
                expected_input_copies,
            "dense T=2/T=4 AVX2 encode retained avoidable input copies");
        require(original == original_before,
            "dense T=2/T=4 AVX2 encode modified caller input");
        if ((test_case.bytes & 63U) == 0 &&
            leo_encode_work_count(test_case.k, test_case.r) != 0)
        {
            require(actual == encode_legacy(original, test_case.r),
                "dense T=2/T=4 AVX2 encode changed legacy parity bytes");
        }
        require(actual == encode_unaligned(
                    avx2_codec.get(), original, test_case.r),
            "dense T=2/T=4 AVX2 encode mishandled unaligned buffers");

        leopard::ff8::TestOnlyResetHighEncodeCounts();
        require_sparse_matches_full(
            encode(avx2_codec.get(), original, test_case.r, true), actual);
        const uint64_t sparse_calls = leopard::ff8::
            TestOnlyGetHighEncodeCounts().small_transform_calls;
        require((sparse_calls != 0) == (test_case.r == 2),
            "T=2/T=4 coarse kernel mishandled a prefix/holey output mask");

        if (test_case.k == 64 && test_case.r == 3)
        {
            leopard::ff8::TestOnlyResetHighEncodeCounts();
            const Shards prefix = encode_prefix(
                avx2_codec.get(), original, test_case.r, 2);
            require(prefix[0] == actual[0] && prefix[1] == actual[1],
                "T=4 coarse kernel changed a requested parity prefix");
            require(leopard::ff8::TestOnlyGetHighEncodeCounts().
                        small_transform_calls == 1,
                "T=4 parity prefix missed the coarse kernel");
        }

        if (scalar.result() == LEO2_SUCCESS)
        {
            Codec scalar_codec(scalar.get(), test_case.k, test_case.r,
                LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8);
            require(actual == encode(
                        scalar_codec.get(), original, test_case.r, false),
                "dense T=2/T=4 AVX2 encode differs from scalar");
        }
        if (ssse3.result() == LEO2_SUCCESS)
        {
            Codec ssse3_codec(ssse3.get(), test_case.k, test_case.r,
                LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8);
            require(actual == encode(
                        ssse3_codec.get(), original, test_case.r, false),
                "dense T=2/T=4 AVX2 encode differs from SSSE3");
        }
        if (avx512.result() == LEO2_SUCCESS)
        {
            Codec avx512_codec(avx512.get(), test_case.k, test_case.r,
                LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8);
            require(actual == encode(
                        avx512_codec.get(), original, test_case.r, false),
                "dense T=2/T=4 AVX2 encode differs from AVX-512");
        }
    }

    // Lock every immediate lower K/byte boundary in the measured staircase.
    // The K<side and T=8 controls also protect the callback's preconditions.
    const TestCase controls[] = {
        { 1, 2, 65536 },
        { 2, 2, 2047 },
        { 3, 4, 2047 },
        { 4, 4, 2047 },
        { 5, 4, 2047 },
        { 6, 3, 2047 },
        { 7, 4, 2047 },
        { 8, 4, 2047 },
        { 9, 4, 2047 },
        { 11, 4, 2047 },
        { 12, 4, 2047 },
        { 16, 4, 2047 },
        { 64, 5, 65536 }
    };
    for (size_t case_i = 0;
         case_i < sizeof(controls) / sizeof(controls[0]); ++case_i)
    {
        const TestCase& test_case = controls[case_i];
        const Shards original = make_original(test_case.k, test_case.bytes);
        Codec codec(avx2.get(), test_case.k, test_case.r,
            LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8);
        leopard::ff8::TestOnlyResetHighEncodeCounts();
        (void)encode(codec.get(), original, test_case.r, false);
        require(leopard::ff8::TestOnlyGetHighEncodeCounts().
                    small_transform_calls == 0,
            "T=2/T=4 coarse-kernel policy escaped its shape bounds");
    }
}

void test_selection_and_bytes(
    Context& automatic,
    Context& scalar,
    Context& ssse3,
    Context& avx2,
    Context& avx512)
{
    Codec auto_codec(automatic.get(), 1000, 200);
    Codec avx2_codec(avx2.get(), 1000, 200);
    Codec avx512_codec(avx512.get(), 1000, 200);

    require(selected_backend(auto_codec.get(), 32, 200, 200) ==
            LEO2_BACKEND_AVX2,
        "AUTO widened a partial transform tile");
    require(selected_backend(auto_codec.get(), 63, 200, 200) ==
            LEO2_BACKEND_AVX2,
        "AUTO widened immediately below the minimum shard length");
    require(selected_backend(auto_codec.get(), 64, 200, 200) ==
            LEO2_BACKEND_AVX512,
        "AUTO did not widen at the complete-tile boundary");
    require(selected_backend(auto_codec.get(), 65, 200, 200) ==
            LEO2_BACKEND_AVX2,
        "AUTO widened an immediate ragged-tail boundary");
    require(selected_backend(auto_codec.get(), 4098, 200, 200) ==
            LEO2_BACKEND_AVX2,
        "AUTO widened an uncalibrated tail");
    require(selected_backend(auto_codec.get(), 4U * 1024U * 1024U,
                200, 200) ==
            LEO2_BACKEND_AVX512,
        "AUTO did not widen a calibrated large complete tile");
    require(selected_backend(auto_codec.get(),
                4U * 1024U * 1024U + 64U, 200, 200) ==
            LEO2_BACKEND_AVX2,
        "AUTO widened immediately above the calibrated shard-length range");
    require(selected_backend(auto_codec.get(), 4096, 199, 200) ==
            LEO2_BACKEND_AVX2,
        "AUTO widened a partial-output encode");
    require(selected_backend(avx2_codec.get(), 4096, 200, 200) ==
            LEO2_BACKEND_AVX2,
        "explicit AVX2 widened");
    require(selected_backend(avx512_codec.get(), 4032, 200, 200) ==
            LEO2_BACKEND_AVX512,
        "explicit AVX512 did not remain exact");
    require(selected_backend(avx512_codec.get(),
                4U * 1024U * 1024U + 64U, 200, 200) ==
            LEO2_BACKEND_AVX512,
        "explicit AVX512 was constrained by AUTO-only byte bounds");

    Codec large(automatic.get(), 4096, 512);
    require(selected_backend(large.get(), 4032, 512, 512) ==
            LEO2_BACKEND_AVX512,
        "large AUTO did not widen a complete tile");
    require(selected_backend(large.get(), 4096, 512, 512) ==
            LEO2_BACKEND_AVX512,
        "large AUTO did not widen in the calibrated cell");
    require(selected_backend(large.get(), 4160, 512, 512) ==
            LEO2_BACKEND_AVX512,
        "large AUTO did not widen a neighboring complete tile");

    Codec minimum_shape(automatic.get(), 8, 2);
    require(leo2_codec_parent_count(minimum_shape.get()) == 16,
        "minimum calibrated parent changed");
    require(selected_backend(minimum_shape.get(), 64, 2, 2) ==
            LEO2_BACKEND_AVX512,
        "minimum K/N transform shape did not widen");
    Codec below_k(automatic.get(), 7, 2);
    require(leo2_codec_parent_count(below_k.get()) == 16,
        "below-K parent does not isolate the K boundary");
    require(selected_backend(below_k.get(), 64, 2, 2) ==
            LEO2_BACKEND_AVX2,
        "AUTO widened immediately below K=8");
    Codec below_parent(automatic.get(), 7, 1);
    require(leo2_codec_parent_count(below_parent.get()) == 8,
        "immediate lower parent is not N=8");
    require(selected_backend(below_parent.get(), 64, 1, 1) ==
            LEO2_BACKEND_AVX2,
        "AUTO widened below N=16");
    Codec direct_xor(automatic.get(), 8, 1);
    require(leo2_codec_parent_count(direct_xor.get()) == 16,
        "R=1 boundary did not retain N=16");
    require(selected_backend(direct_xor.get(), 4096, 1, 1) ==
            LEO2_BACKEND_AVX2,
        "T=1 direct codec unnecessarily qualified transform widening");

    Codec maximum_r(automatic.get(), 8, 4096);
    require(selected_backend(maximum_r.get(), 64, 4096, 4096) ==
            LEO2_BACKEND_AVX512,
        "AUTO did not widen at R=4096");
    Codec above_r(automatic.get(), 8, 4097);
    require(selected_backend(above_r.get(), 64, 4097, 4097) ==
            LEO2_BACKEND_AVX2,
        "AUTO widened immediately above R=4096");

    Codec low_profile(automatic.get(), 8, 8, LEO2_PROFILE_LOW_V1);
    require(selected_backend(low_profile.get(), 64, 8, 8) ==
            LEO2_BACKEND_AVX2,
        "AUTO widened a non-legacy-high profile");
    Codec gf8_unbalanced(automatic.get(), 8, 2,
        LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8);
    require(selected_backend(gf8_unbalanced.get(), 4096, 2, 2) ==
            LEO2_BACKEND_AVX2,
        "AUTO widened an unqualified GF8 codec");

    static const uint32_t gf8_sides[] = { 8, 16, 32, 64 };
    for (size_t side_i = 0;
         side_i < sizeof(gf8_sides) / sizeof(gf8_sides[0]); ++side_i)
    {
        const uint32_t side = gf8_sides[side_i];
        Codec balanced(automatic.get(), side, side,
            LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8);
        const size_t execution_bytes = side == 8 ? 8192 : 4096;
        if (side == 8)
        {
            require(selected_backend(balanced.get(), 1984, side, side) ==
                    LEO2_BACKEND_AVX2,
                "balanced GF8 T=8 widened below the 2 KiB singleton");
            require(selected_backend(balanced.get(), 2048, side, side) ==
                    LEO2_BACKEND_AVX512,
                "balanced GF8 T=8 did not widen at 2 KiB");
            require(selected_backend(balanced.get(), 2112, side, side) ==
                    LEO2_BACKEND_AVX2,
                "balanced GF8 T=8 widened above the 2 KiB singleton");
            require(selected_backend(balanced.get(), 4096, side, side) ==
                    LEO2_BACKEND_AVX2,
                "balanced GF8 T=8 ignored the 4 KiB exact-main veto");
            require(selected_backend(balanced.get(), 8128, side, side) ==
                    LEO2_BACKEND_AVX2,
                "balanced GF8 T=8 widened below the upper interval");
            require(selected_backend(balanced.get(), 8192, side, side) ==
                    LEO2_BACKEND_AVX512,
                "balanced GF8 T=8 did not widen at 8 KiB");
        }
        else
        {
            const uint64_t minimum_bytes = side == 16 ? 4096 : 2048;
            require(selected_backend(
                        balanced.get(), minimum_bytes - 64, side, side) ==
                    LEO2_BACKEND_AVX2,
                "balanced GF8 AUTO widened below the calibrated byte range");
            require(selected_backend(
                        balanced.get(), minimum_bytes, side, side) ==
                    LEO2_BACKEND_AVX512,
                "balanced GF8 AUTO did not widen at the lower byte boundary");
        }
        if (side == 16)
            require(selected_backend(balanced.get(), 2048, side, side) ==
                    LEO2_BACKEND_AVX2,
                "balanced GF8 T=16 widened at the inconclusive 2 KiB cell");
        require(selected_backend(
                    balanced.get(), execution_bytes, side, side) ==
                LEO2_BACKEND_AVX512,
            "balanced GF8 AUTO did not widen in a calibrated cell");
        require(selected_backend(balanced.get(), 65536, side, side) ==
                LEO2_BACKEND_AVX512,
            "balanced GF8 AUTO did not widen at the upper byte boundary");
        require(selected_backend(balanced.get(), 65600, side, side) ==
                LEO2_BACKEND_AVX2,
            "balanced GF8 AUTO widened above the calibrated byte range");
        require(selected_backend(balanced.get(), 4097, side, side) ==
                LEO2_BACKEND_AVX2,
            "balanced GF8 AUTO widened a ragged tail");
        require(selected_backend(
                    balanced.get(), execution_bytes, side - 1, side) ==
                LEO2_BACKEND_AVX2,
            "balanced GF8 AUTO widened a partial-output encode");

        Codec explicit_avx2(avx2.get(), side, side,
            LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8);
        Codec explicit_avx512(avx512.get(), side, side,
            LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8);
        Codec explicit_scalar(scalar.get(), side, side,
            LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8);
        Codec explicit_ssse3(ssse3.get(), side, side,
            LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8);
        require(selected_backend(explicit_avx2.get(), 4096, side, side) ==
                LEO2_BACKEND_AVX2,
            "balanced GF8 explicit AVX2 widened");
        require(selected_backend(explicit_avx512.get(), 4097, side, side) ==
                LEO2_BACKEND_AVX512,
            "balanced GF8 explicit AVX512 was constrained by AUTO bounds");

        Shards balanced_original = make_original(side, execution_bytes);
        leopard::ff8::TestOnlyResetHighEncodeCounts();
        const Shards balanced_auto = encode(
            balanced.get(), balanced_original, side, false);
        const leopard::ff8::TestOnlyHighEncodeCounts whole_counts =
            leopard::ff8::TestOnlyGetHighEncodeCounts();
        require(whole_counts.whole_transform_calls == 1,
            "balanced GF8 AUTO did not execute the coarse transform");
        const Shards balanced_avx2 = encode(
            explicit_avx2.get(), balanced_original, side, false);
        const Shards balanced_avx512 = encode(
            explicit_avx512.get(), balanced_original, side, false);
        require(balanced_auto == balanced_avx2 &&
                balanced_auto == balanced_avx512 &&
                balanced_auto == encode(
                    explicit_scalar.get(), balanced_original, side, false) &&
                balanced_auto == encode(
                    explicit_ssse3.get(), balanced_original, side, false) &&
                balanced_auto == encode_legacy(balanced_original, side),
            "balanced GF8 coarse transform changed legacy parity bytes");
        if (side == 8)
            require_encode_overlap_rejected(
                balanced.get(), balanced_original, side);
        leopard::ff8::TestOnlyResetHighEncodeCounts();
        require_sparse_matches_full(
            encode(balanced.get(), balanced_original, side, true),
            balanced_auto);
        require(leopard::ff8::TestOnlyGetHighEncodeCounts().
                    whole_transform_calls == 0,
            "balanced GF8 sparse encode used the dense coarse transform");

        static const size_t tail_bytes[] = { 1025, 2049, 4097 };
        for (size_t tail_i = 0;
             tail_i < sizeof(tail_bytes) / sizeof(tail_bytes[0]); ++tail_i)
        {
            const Shards tail_original =
                make_original(side, tail_bytes[tail_i]);
            const Shards tail_reference = encode(
                explicit_scalar.get(), tail_original, side, false);
            const Shards tail_ssse3 = encode(
                explicit_ssse3.get(), tail_original, side, false);
            const Shards tail_avx2 = encode(
                explicit_avx2.get(), tail_original, side, false);
            leopard::ff8::TestOnlyResetHighEncodeCounts();
            const Shards tail_avx512 = encode(
                explicit_avx512.get(), tail_original, side, false);
            const uint64_t expected_whole_calls = tail_bytes[tail_i] >= 2049
                ? 1 : 0;
            require(leopard::ff8::TestOnlyGetHighEncodeCounts().
                    whole_transform_calls == expected_whole_calls,
                "balanced GF8 tail used the wrong coarse-transform split");
            require(tail_reference == tail_ssse3 &&
                    tail_reference == tail_avx2 &&
                    tail_reference == tail_avx512,
                "balanced GF8 coarse transform changed a byte tail");
            require(tail_reference == encode_unaligned(
                        explicit_avx2.get(), tail_original, side) &&
                    tail_reference == encode_unaligned(
                        explicit_avx512.get(), tail_original, side),
                "balanced GF8 coarse transform mishandled unaligned buffers");
        }
    }

    // Exercise concurrent immutable-codec use through the exact T=8 AUTO
    // route independently of the larger neighboring callback below.
    Codec concurrent_t8(automatic.get(), 8, 8,
        LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8);
    const Shards concurrent_t8_original = make_original(8, 8192);
    const Shards concurrent_t8_reference = encode(
        concurrent_t8.get(), concurrent_t8_original, 8, false);
    std::atomic<unsigned> concurrent_t8_failures(0);
    std::vector<std::thread> concurrent_t8_threads;
    for (unsigned lane = 0; lane < 4; ++lane)
    {
        concurrent_t8_threads.push_back(std::thread([&]() {
            try
            {
                if (encode(concurrent_t8.get(), concurrent_t8_original,
                        8, false) != concurrent_t8_reference)
                    concurrent_t8_failures.fetch_add(
                        1, std::memory_order_relaxed);
            }
            catch (...)
            {
                concurrent_t8_failures.fetch_add(
                    1, std::memory_order_relaxed);
            }
        }));
    }
    for (size_t i = 0; i < concurrent_t8_threads.size(); ++i)
        concurrent_t8_threads[i].join();
    require(concurrent_t8_failures.load(std::memory_order_relaxed) == 0,
        "concurrent exact T=8 AUTO encode failed");

    Codec gf8_t8(automatic.get(), 8, 8,
        LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8);
    Codec gf8_above_side(automatic.get(), 128, 128,
        LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8);
    Codec gf8_shortened(automatic.get(), 15, 16,
        LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8);
    Codec gf8_punctured(automatic.get(), 16, 15,
        LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8);
    require(selected_backend(gf8_t8.get(), 4096, 8, 8) ==
            LEO2_BACKEND_AVX2,
        "balanced GF8 AUTO ignored the T=8 4 KiB veto");
    require(selected_backend(gf8_t8.get(), 8192, 8, 8) ==
            LEO2_BACKEND_AVX512,
        "balanced GF8 AUTO did not widen qualified T=8");
    require(selected_backend(gf8_above_side.get(), 4096, 128, 128) ==
            LEO2_BACKEND_AVX2,
        "balanced GF8 AUTO widened T=128");
    require(selected_backend(gf8_shortened.get(), 4096, 16, 16) ==
            LEO2_BACKEND_AVX2,
        "balanced GF8 T=16 widened an unqualified shortened input block");
    require(selected_backend(gf8_punctured.get(), 4096, 15, 15) ==
            LEO2_BACKEND_AVX2,
        "balanced GF8 T=16 widened an unqualified punctured parity block");

    struct GF8CoarseCase
    {
        uint32_t k;
        uint32_t r;
    };
    static const GF8CoarseCase coarse_cases[] = {
        { 8, 8 }, { 7, 8 }, { 8, 7 }, { 7, 7 },
        { 15, 16 }, { 16, 15 }, { 15, 15 },
        { 31, 32 }, { 32, 31 }, { 31, 31 },
        { 63, 64 }, { 64, 63 }, { 63, 63 }
    };
    static const size_t coarse_bytes[] = {
        2048, 2049, 4096, 4097, 8192, 65536, 65537
    };
    for (size_t case_i = 0;
         case_i < sizeof(coarse_cases) / sizeof(coarse_cases[0]); ++case_i)
    {
        const GF8CoarseCase& current = coarse_cases[case_i];
        Codec automatic_codec(automatic.get(), current.k, current.r,
            LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8);
        Codec scalar_codec(scalar.get(), current.k, current.r,
            LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8);
        Codec ssse3_codec(ssse3.get(), current.k, current.r,
            LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8);
        Codec avx2_codec(avx2.get(), current.k, current.r,
            LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8);
        Codec avx512_codec(avx512.get(), current.k, current.r,
            LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8);

        const uint32_t side = current.r <= 8 ? 8U :
            (current.r <= 16 ? 16U : (current.r <= 32 ? 32U : 64U));
        const bool qualified_shape =
            (side == 8 && current.k == 8 && current.r == 8) ||
            (side == 64 &&
                (current.k == 64 || current.k == 63) &&
                (current.r == 64 || current.r == 63));

        // AUTO widens only the exact aligned cells that passed the isolated
        // crossover gate.  Explicit AVX-512 exercises every neighboring
        // candidate independently of the production promotion decision.
        for (size_t byte_i = 0;
             byte_i < sizeof(coarse_bytes) / sizeof(coarse_bytes[0]); ++byte_i)
        {
            const size_t bytes = coarse_bytes[byte_i];
            const bool qualified_bytes = side == 8
                ? bytes == 2U * 1024U ||
                    (bytes >= 8U * 1024U && bytes <= 64U * 1024U)
                : current.r == 64 || bytes == 64U * 1024U;
            const leo2_backend expected_auto = qualified_shape &&
                    qualified_bytes && (bytes & 63U) == 0
                ? LEO2_BACKEND_AVX512 : LEO2_BACKEND_AVX2;
            require(selected_backend(automatic_codec.get(), bytes,
                        current.r, current.r) == expected_auto,
                "GF8 coarse neighbor selected the wrong AUTO backend");
            require(selected_backend(avx512_codec.get(), bytes,
                        current.r, current.r) == LEO2_BACKEND_AVX512,
                "explicit GF8 coarse neighbor lost AVX-512");

            const Shards original = make_original(current.k, bytes);
            const Shards reference = encode(
                scalar_codec.get(), original, current.r, false);
            require(reference == encode(
                        ssse3_codec.get(), original, current.r, false) &&
                    reference == encode(
                        avx2_codec.get(), original, current.r, false),
                "GF8 coarse neighbor changed a mature backend result");

            leopard::ff8::TestOnlyResetHighEncodeCounts();
            const Shards widened = encode(
                avx512_codec.get(), original, current.r, false);
            require(leopard::ff8::TestOnlyGetHighEncodeCounts().
                        whole_transform_calls == 1,
                "GF8 coarse neighbor did not execute one aligned-prefix callback");
            const Shards automatic_result = encode(
                automatic_codec.get(), original, current.r, false);
            require(reference == widened && reference == automatic_result,
                "GF8 coarse neighbor changed parity bytes");
            const uint64_t actual_whole_transform_calls =
                leopard::ff8::TestOnlyGetHighEncodeCounts().
                    whole_transform_calls;
            /*
                The explicit AVX-512 reference contributes one callback.
                AUTO contributes another either when it selects AVX-512 for
                the aligned body or when exact T=8 AVX2 encodes a staged
                64-byte tail through its independently qualified callback.
            */
            const uint64_t expected_whole_transform_calls =
                1U + (
                    expected_auto == LEO2_BACKEND_AVX512 ||
                    (current.k == 8 && current.r == 8 &&
                     (bytes & 63U) != 0)
                        ? 1U : 0U);
            if (actual_whole_transform_calls !=
                expected_whole_transform_calls)
            {
                throw std::runtime_error(
                    "GF8 coarse neighbor executed the wrong AUTO callback "
                    "route: K=" + std::to_string(current.k) +
                    " R=" + std::to_string(current.r) +
                    " bytes=" + std::to_string(bytes) +
                    " expected=" +
                    std::to_string(expected_whole_transform_calls) +
                    " actual=" +
                    std::to_string(actual_whole_transform_calls));
            }
            if ((bytes & 63U) == 0 && current.r <= current.k)
                require(reference == encode_legacy(original, current.r),
                    "GF8 coarse neighbor differs from legacy Leopard");

            if ((bytes & 63U) != 0)
            {
                require(reference == encode_unaligned(
                            avx2_codec.get(), original, current.r) &&
                        reference == encode_unaligned(
                            avx512_codec.get(), original, current.r),
                    "GF8 coarse neighbor mishandled an unaligned byte tail");
            }
        }

        const Shards partial_original = make_original(current.k, 4097);
        const Shards partial_reference = encode(
            scalar_codec.get(), partial_original, current.r, false);
        leopard::ff8::TestOnlyResetHighEncodeCounts();
        require_sparse_matches_full(encode(
                avx512_codec.get(), partial_original, current.r, true),
            partial_reference);
        require(leopard::ff8::TestOnlyGetHighEncodeCounts().
                    whole_transform_calls == 0,
            "partial-output GF8 neighbor used the dense coarse callback");
    }

    // One ragged, shortened-and-punctured codec is also executed concurrently;
    // each call owns its output and scratch while sharing the immutable codec.
    Codec concurrent_codec(avx512.get(), 63, 63,
        LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8);
    const Shards concurrent_original = make_original(63, 4097);
    const Shards concurrent_reference = encode(
        concurrent_codec.get(), concurrent_original, 63, false);
    std::atomic<unsigned> coarse_failures(0);
    std::vector<std::thread> coarse_threads;
    for (unsigned lane = 0; lane < 4; ++lane)
    {
        coarse_threads.push_back(std::thread([&]() {
            try
            {
                if (encode(concurrent_codec.get(), concurrent_original,
                        63, false) != concurrent_reference)
                    coarse_failures.fetch_add(1, std::memory_order_relaxed);
            }
            catch (...)
            {
                coarse_failures.fetch_add(1, std::memory_order_relaxed);
            }
        }));
    }
    for (size_t i = 0; i < coarse_threads.size(); ++i)
        coarse_threads[i].join();
    require(coarse_failures.load(std::memory_order_relaxed) == 0,
        "concurrent GF8 coarse neighbor encode failed");

    Codec maximum_parent(automatic.get(), 60000, 1024);
    require(selected_backend(maximum_parent.get(), 1024, 1024, 1024) ==
            LEO2_BACKEND_AVX512,
        "maximum-parent GF16 high transform did not widen");

    const Shards original = make_original(1000, 4096);
    const Shards automatic_parity = encode(
        auto_codec.get(), original, 200, false);
    const Shards avx2_parity = encode(
        avx2_codec.get(), original, 200, false);
    const Shards avx512_parity = encode(
        avx512_codec.get(), original, 200, false);
    require(automatic_parity == avx2_parity &&
            automatic_parity == avx512_parity,
        "AUTO widening changed GF16 legacy parity bytes");
    require_sparse_matches_full(
        encode(auto_codec.get(), original, 200, true), automatic_parity);

    std::atomic<unsigned> failures(0);
    std::vector<std::thread> threads;
    for (unsigned lane = 0; lane < 4; ++lane)
    {
        threads.push_back(std::thread([&]() {
            try
            {
                require(encode(auto_codec.get(), original, 200, false) ==
                        automatic_parity,
                    "concurrent AUTO parity mismatch");
            }
            catch (...)
            {
                failures.fetch_add(1, std::memory_order_relaxed);
            }
        }));
    }
    for (size_t i = 0; i < threads.size(); ++i)
        threads[i].join();
    require(failures.load(std::memory_order_relaxed) == 0,
        "concurrent AUTO widening failed");
}

} // namespace

int main()
{
    try
    {
        require(leo_init() == Leopard_Success, "Leopard initialization");
        Context automatic(LEO2_BACKEND_AUTO);
        require_result(automatic.result(), LEO2_SUCCESS, "AUTO context");
        Context scalar(LEO2_BACKEND_SCALAR);
        Context ssse3(LEO2_BACKEND_SSSE3);
        Context avx2(LEO2_BACKEND_AVX2);
        Context avx512(LEO2_BACKEND_AVX512);
        // Explicit-only, like AVX-512: an unqualified host returns
        // LEO2_UNSUPPORTED and require_explicit_backend returns early.
        Context gfni(LEO2_BACKEND_GFNI);

        if (std::getenv("LEO2_TEST_T8_ONE_BLOCK_ONLY"))
        {
            test_t8_one_block_mixed_binding(avx2);
            test_t8_two_block_mixed_binding(avx2);
            std::puts("Leopard2 T8 mixed binding passed");
            return 0;
        }

        require_explicit_backend(scalar, LEO2_BACKEND_SCALAR);
        require_explicit_backend(ssse3, LEO2_BACKEND_SSSE3);
        require_explicit_backend(avx2, LEO2_BACKEND_AVX2);
        require_explicit_backend(avx512, LEO2_BACKEND_AVX512);
        require_explicit_backend(gfni, LEO2_BACKEND_GFNI);
        const uint64_t detected_l3_bytes =
            leopard::backend::DetectConservativeL3Bytes();
        require_context_cache_policy_wiring(
            automatic, detected_l3_bytes, false);
        require_context_cache_policy_wiring(
            scalar, detected_l3_bytes, false);
        require_context_cache_policy_wiring(
            ssse3, detected_l3_bytes, false);
        require_context_cache_policy_wiring(
            avx2, detected_l3_bytes, true);
        require_context_cache_policy_wiring(
            avx512, detected_l3_bytes, false);
        require_context_cache_policy_wiring(
            gfni, detected_l3_bytes, false);
        test_balanced_execution_tile_geometry();
        test_linux_cache_topology();
        test_gf16_cache_policy(avx2);
        test_gf16_cache_policy_bytes(avx2);
        test_t8_one_block_mixed_binding(avx2);
        test_t8_two_block_mixed_binding(avx2);

        test_small_high_encode(scalar, ssse3, avx2, avx512);

        if (leopard::backend::IsCalibratedAutoAVX512EncodeHost() &&
            leo2_context_backend(automatic.get()) == LEO2_BACKEND_AVX2 &&
            scalar.result() == LEO2_SUCCESS &&
            ssse3.result() == LEO2_SUCCESS &&
            avx2.result() == LEO2_SUCCESS &&
            avx512.result() == LEO2_SUCCESS)
        {
            test_selection_and_bytes(
                automatic, scalar, ssse3, avx2, avx512);
        }
        else
        {
            Codec codec(automatic.get(), 1000, 200);
            require(selected_backend(codec.get(), 4096, 200, 200) ==
                    leo2_context_backend(automatic.get()),
                "unsupported host widened away from its AUTO baseline");
        }

        std::printf("Leopard2 AUTO encode backend passed\n");
        return 0;
    }
    catch (const std::exception& error)
    {
        std::fprintf(stderr, "Leopard2 AUTO encode backend failed: %s\n",
            error.what());
        return 1;
    }
}
