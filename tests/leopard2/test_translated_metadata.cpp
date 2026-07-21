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

#include "Leopard2Direct.h"
#include "Leopard2Dispatch.h"

#include <cstring>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#ifndef LEO2_EXPECT_GF8
#define LEO2_EXPECT_GF8 1
#endif

#ifndef LEO2_EXPECT_GF16
#define LEO2_EXPECT_GF16 1
#endif

#ifndef LEO2_EXPECT_HOOK_NATIVE_METADATA
#define LEO2_EXPECT_HOOK_NATIVE_METADATA 0
#endif

namespace {

void require(bool condition, const std::string& message)
{
    if (!condition)
        throw std::runtime_error(message);
}

void require_success(leo2_result result, const std::string& operation)
{
    if (result == LEO2_SUCCESS)
        return;
    throw std::runtime_error(operation + ": " + leo2_result_string(result));
}

struct Shape
{
    leo2_field field;
    uint32_t k;
    uint32_t r;
    uint32_t side;
    uint32_t parent;
    size_t element_bytes;
    const char* name;
};

leo2_codec* create_codec(
    leo2_context* context,
    const Shape& shape,
    leo2_profile profile,
    uint32_t flags)
{
    leo2_codec_options options;
    memset(&options, 0, sizeof(options));
    options.struct_size = sizeof(options);
    options.flags = flags;
    leo2_codec* codec = NULL;
    require_success(leo2_codec_create(context, shape.k, shape.r,
        profile, shape.field, &options, &codec), "codec create");
    require(codec != NULL, "codec create returned null");
    return codec;
}

leopard2_internal::CodecDecodeMetadataInfo metadata(const leo2_codec* codec)
{
    leopard2_internal::CodecDecodeMetadataInfo info;
    require(leopard2_internal::GetCodecDecodeMetadataInfo(codec, &info),
        "metadata introspection failed");
    return info;
}

void require_empty_translated(
    const leopard2_internal::CodecDecodeMetadataInfo& info,
    const std::string& label)
{
    require(info.translated_permanent_erased_bytes == 0 &&
            info.translated_locator_bytes == 0 &&
            info.translated_factor_bytes == 0,
        label + " retained translated metadata");
}

void test_translated_codec_metadata(
    leo2_context* context,
    const Shape& shape,
    uint32_t flags)
{
    leo2_codec* codec = create_codec(context, shape,
        LEO2_PROFILE_LEGACY_HIGH_V1, flags);
    const leopard2_internal::CodecDecodeMetadataInfo info = metadata(codec);
    const std::string label = std::string(shape.name) + " translated codec";
    require(info.permanent_erased_bytes == shape.parent,
        label + " lost its parent erasure mask");
    require(info.translated_permanent_erased_bytes == shape.parent,
        label + " did not retain its translated erasure mask");
    require(info.translated_locator_bytes ==
            shape.parent * shape.element_bytes,
        label + " did not retain its reachable translated locator");
    require(info.translated_factor_bytes == shape.element_bytes,
        label + " did not retain its one reachable Algorithm-4 factor");
#if LEO2_EXPECT_HOOK_NATIVE_METADATA
    require(info.native_locator_bytes == shape.parent * shape.element_bytes,
        label + " hook build lost its native locator control");
    require(info.native_factor_bytes == shape.parent * shape.element_bytes,
        label + " hook build lost its Algorithm-5 factors");
    require_success(leo2_test_codec_set_decode_mode(codec,
        LEO2_TEST_DECODE_FORCE_NATIVE_HIGH),
        label + " select native-high hook");
#else
    require(info.native_locator_bytes == 0,
        label + " retained unreachable native locator storage");
    require(info.native_factor_bytes == 0,
        label + " retained unreachable Algorithm-5 factors");
#endif
    if (flags == 0)
    {
        std::cout << "{\"build\":\""
#if LEO2_EXPECT_HOOK_NATIVE_METADATA
                  << "hooks"
#else
                  << "production"
#endif
                  << "\",\"field\":\"" << shape.name
                  << "\",\"parent\":" << shape.parent
                  << ",\"permanent_erased_bytes\":"
                  << info.permanent_erased_bytes
                  << ",\"native_locator_bytes\":"
                  << info.native_locator_bytes
                  << ",\"native_factor_bytes\":"
                  << info.native_factor_bytes
                  << ",\"translated_permanent_erased_bytes\":"
                  << info.translated_permanent_erased_bytes
                  << ",\"translated_locator_bytes\":"
                  << info.translated_locator_bytes
                  << ",\"translated_factor_bytes\":"
                  << info.translated_factor_bytes << "}\n";
    }
    leo2_codec_destroy(codec);
}

void test_forced_generic_metadata(leo2_context* context, const Shape& shape)
{
    leo2_codec* codec = create_codec(context, shape,
        LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_CODEC_FORCE_GENERIC_DECODE);
    const leopard2_internal::CodecDecodeMetadataInfo info = metadata(codec);
    const std::string label = std::string(shape.name) + " forced generic";
    require(info.permanent_erased_bytes == shape.parent,
        label + " lost its parent erasure mask");
    require(info.native_locator_bytes == shape.parent * shape.element_bytes,
        label + " lost its reachable native locator");
    require(info.native_factor_bytes == 0,
        label + " prepared specialized factors");
    require_empty_translated(info, label);
    leo2_codec_destroy(codec);
}

void test_nontranslated_profiles(leo2_context* context, const Shape& shape)
{
    Shape high = shape;
    high.k = shape.side + 1;
    high.r = shape.r;
    high.parent = shape.parent * 2;
    leo2_codec* codec = create_codec(context, high,
        LEO2_PROFILE_LEGACY_HIGH_V1, 0);
    leopard2_internal::CodecDecodeMetadataInfo info = metadata(codec);
    const std::string high_label = std::string(shape.name) +
        " nontranslated high";
    require(info.native_factor_bytes == high.parent * shape.element_bytes,
        high_label + " lost reachable Algorithm-5 factors");
    require_empty_translated(info, high_label);
    leo2_codec_destroy(codec);

    codec = create_codec(context, shape, LEO2_PROFILE_LOW_V1, 0);
    info = metadata(codec);
    const std::string low_label = std::string(shape.name) + " low profile";
    require(info.native_factor_bytes == shape.element_bytes,
        low_label + " lost its reachable Algorithm-4 factor");
    require(info.native_locator_bytes == shape.parent * shape.element_bytes,
        low_label + " lost its reachable locator");
    require_empty_translated(info, low_label);
    leo2_codec_destroy(codec);
}

leopard2_internal::DecodePathInfo plan_path(
    const leo2_codec* codec,
    uint32_t losses)
{
    const uint32_t k = leo2_codec_original_count(codec);
    const uint32_t r = leo2_codec_recovery_count(codec);
    std::vector<uint8_t> original_present(k, 1);
    std::vector<uint8_t> recovery_present(r, 1);
    for (uint32_t i = 0; i < losses; ++i)
        original_present[i] = 0;
    leo2_decode_plan* plan = NULL;
    require_success(leo2_decode_plan_create(codec, original_present.data(),
        recovery_present.data(), &plan), "plan create");
    leopard2_internal::DecodePathInfo path;
    require_success(leopard2_internal::GetDecodePlanPathInfo(
        plan, 64, false, &path), "plan path");
    leo2_decode_plan_destroy(plan);
    return path;
}

void test_plan_reachability(leo2_context* context, const Shape& shape)
{
    leo2_codec* automatic = create_codec(context, shape,
        LEO2_PROFILE_LEGACY_HIGH_V1, 0);
    require(plan_path(automatic, 0).path ==
            leopard2_internal::kDecodePathNoOp,
        std::string(shape.name) + " no-loss plan was not a no-op");
#if LEO2_EXPECT_HOOK_NATIVE_METADATA
    require_success(leo2_test_codec_set_decode_mode(automatic,
        LEO2_TEST_DECODE_AUTO), "select hook AUTO");
#endif
    require(plan_path(automatic, 1).path ==
            leopard2_internal::kDecodePathDirect,
        std::string(shape.name) + " AUTO one-loss plan missed direct repair");
    const leopard2_internal::DecodePathInfo translated =
        plan_path(automatic, 5);
    require(translated.rule == leopard2_internal::kDecodeRuleTranslatedLow,
        std::string(shape.name) +
            " AUTO transform plan did not capture translated Algorithm 4");
#if LEO2_EXPECT_HOOK_NATIVE_METADATA
    require_success(leo2_test_codec_set_decode_mode(automatic,
        LEO2_TEST_DECODE_FORCE_NATIVE_HIGH), "select native-high hook");
    const leopard2_internal::DecodePathInfo native = plan_path(automatic, 1);
    require(native.path == leopard2_internal::kDecodePathMaterialized ||
            native.path == leopard2_internal::kDecodePathTiled,
        std::string(shape.name) + " native-high hook missed transform path");
    require(native.rule != leopard2_internal::kDecodeRuleTranslatedLow,
        std::string(shape.name) + " native-high hook remained translated");
#endif
    leo2_codec_destroy(automatic);

    const uint32_t translated_flags[] = {
        LEO2_CODEC_FORCE_SPECIALIZED_DECODE,
        LEO2_CODEC_FORCE_TILED_DECODE,
        LEO2_CODEC_FORCE_MATERIALIZED_DECODE
    };
    for (size_t i = 0;
         i < sizeof(translated_flags) / sizeof(translated_flags[0]); ++i)
    {
        leo2_codec* forced = create_codec(context, shape,
            LEO2_PROFILE_LEGACY_HIGH_V1, translated_flags[i]);
        require(plan_path(forced, 1).rule ==
                leopard2_internal::kDecodeRuleTranslatedLow,
            std::string(shape.name) +
                " forced specialized/workspace plan was not translated");
        leo2_codec_destroy(forced);
    }

    leo2_codec* generic = create_codec(context, shape,
        LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_CODEC_FORCE_GENERIC_DECODE);
    const leopard2_internal::DecodePathInfo generic_path =
        plan_path(generic, 1);
    require(generic_path.path == leopard2_internal::kDecodePathGeneric &&
            generic_path.rule == leopard2_internal::kDecodeRuleForcedGeneric,
        std::string(shape.name) + " forced generic plan was not generic");
    leo2_codec_destroy(generic);
}

void run_shape(leo2_context* context, const Shape& shape)
{
    const uint32_t flags[] = {
        0,
        LEO2_CODEC_FORCE_SPECIALIZED_DECODE,
        LEO2_CODEC_FORCE_TILED_DECODE,
        LEO2_CODEC_FORCE_MATERIALIZED_DECODE
    };
    for (size_t i = 0; i < sizeof(flags) / sizeof(flags[0]); ++i)
        test_translated_codec_metadata(context, shape, flags[i]);
    test_forced_generic_metadata(context, shape);
    test_nontranslated_profiles(context, shape);
    test_plan_reachability(context, shape);
}

} // namespace

int main()
{
    try
    {
        leo2_context_options options;
        memset(&options, 0, sizeof(options));
        options.struct_size = sizeof(options);
        options.backend = LEO2_BACKEND_SCALAR;
        options.thread_count = 1;
        leo2_context* context = NULL;
        require_success(leo2_context_create(&options, &context),
            "context create");
#if LEO2_EXPECT_GF8
        const Shape gf8 = {
            LEO2_FIELD_GF8, 9, 9, 16, 32, 1, "GF8"
        };
        run_shape(context, gf8);
#endif
#if LEO2_EXPECT_GF16
        const Shape gf16 = {
            LEO2_FIELD_GF16, 9, 9, 16, 32, 2, "GF16"
        };
        run_shape(context, gf16);
#endif
        require(!leopard2_internal::GetCodecDecodeMetadataInfo(NULL, NULL),
            "null metadata query succeeded");
        leo2_context_destroy(context);
        std::cout << "translated metadata structural tests passed\n";
        return 0;
    }
    catch (const std::exception& error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
