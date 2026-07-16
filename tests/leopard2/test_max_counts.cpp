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

#include "leopard2.h"

#include <array>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <new>
#include <stdexcept>
#include <string>
#include <vector>

#if defined(_MSC_VER)
#include <malloc.h>
#endif

namespace {

typedef std::array<unsigned char, 2> Word;

class AlignedBuffer
{
public:
    explicit AlignedBuffer(size_t bytes)
        : data_(NULL)
        , bytes_(bytes)
    {
        if (bytes == 0)
            return;
#if defined(_MSC_VER)
        data_ = _aligned_malloc(bytes, leo2_scratch_alignment());
#else
        if (posix_memalign(&data_, leo2_scratch_alignment(), bytes) != 0)
            data_ = NULL;
#endif
        if (!data_)
            throw std::bad_alloc();
        memset(data_, 0, bytes);
    }

    ~AlignedBuffer()
    {
#if defined(_MSC_VER)
        _aligned_free(data_);
#else
        free(data_);
#endif
    }

    AlignedBuffer(AlignedBuffer&& other)
        : data_(other.data_)
        , bytes_(other.bytes_)
    {
        other.data_ = NULL;
        other.bytes_ = 0;
    }

    void* data() { return data_; }
    size_t size() const { return bytes_; }

private:
    AlignedBuffer(const AlignedBuffer&);
    AlignedBuffer& operator=(const AlignedBuffer&);

    void* data_;
    size_t bytes_;
};

void require(bool condition, const char* operation)
{
    if (!condition)
        throw std::runtime_error(operation);
}

void require_result(leo2_result result, const char* operation)
{
    if (result != LEO2_SUCCESS)
    {
        throw std::runtime_error(
            std::string(operation) + ": " + leo2_result_string(result));
    }
}

leo2_codec* make_codec(
    leo2_context* context,
    uint32_t k,
    uint32_t r,
    leo2_profile profile)
{
    leo2_codec* codec = NULL;
    require_result(leo2_codec_create(context, k, r, profile,
        LEO2_FIELD_GF16, NULL, &codec), "maximum codec create");
    require(codec != NULL, "maximum codec create returned null");
    require(leo2_codec_parent_count(codec) == 65536,
        "maximum codec parent is not the full GF16 field");
    return codec;
}

AlignedBuffer encode_scratch(const leo2_codec* codec)
{
    size_t bytes = 0;
    require_result(leo2_encode_scratch_size(codec, 2, &bytes),
        "maximum encode scratch query");
    return AlignedBuffer(bytes);
}

void recover_missing(
    const leo2_codec* codec,
    const std::vector<uint8_t>& original_present,
    const std::vector<uint8_t>& recovery_present,
    const std::vector<const void*>& original,
    const std::vector<const void*>& recovery,
    std::vector<void*>& restored)
{
    leo2_decode_plan* plan = NULL;
    require_result(leo2_decode_plan_create(codec, &original_present[0],
        &recovery_present[0], &plan), "maximum decode plan create");
    size_t scratch_bytes = 0;
    require_result(leo2_decode_plan_scratch_size(plan, 2, &scratch_bytes),
        "maximum decode scratch query");
    AlignedBuffer scratch(scratch_bytes);
    require_result(leo2_decode_plan_execute(plan, 2, &original[0],
        &recovery[0], &restored[0], scratch.data(), scratch.size()),
        "maximum decode execute");
    leo2_decode_plan_destroy(plan);
}

void test_high_r1_limit(leo2_context* context, uint64_t& compared_bytes)
{
    const uint32_t k = 65535;
    leo2_codec* codec = make_codec(
        context, k, 1, LEO2_PROFILE_LEGACY_HIGH_V1);
    require(leo2_codec_padded_side(codec) == 1,
        "maximum high R=1 padded side mismatch");

    std::vector<Word> original(k);
    std::vector<const void*> original_input(k);
    Word expected = {{ 0, 0 }};
    for (uint32_t i = 0; i < k; ++i)
    {
        original[i][0] = static_cast<unsigned char>(i * 17u + 3u);
        original[i][1] = static_cast<unsigned char>(i * 29u + (i >> 8));
        expected[0] ^= original[i][0];
        expected[1] ^= original[i][1];
        original_input[i] = original[i].data();
    }

    Word parity = {{ 0, 0 }};
    void* parity_output[] = { parity.data() };
    AlignedBuffer scratch = encode_scratch(codec);
    require_result(leo2_encode(codec, 2, &original_input[0], parity_output,
        scratch.data(), scratch.size()), "maximum high R=1 encode");
    require(parity == expected, "maximum high R=1 parity is not XOR");
    compared_bytes += 2;

    const uint32_t missing = 32768;
    std::vector<uint8_t> original_present(k, 1);
    original_present[missing] = 0;
    std::vector<uint8_t> recovery_present(1, 1);
    original_input[missing] = NULL;
    std::vector<const void*> recovery_input(1, parity.data());
    Word restored_word = {{ 0, 0 }};
    std::vector<void*> restored(k, NULL);
    restored[missing] = restored_word.data();
    recover_missing(codec, original_present, recovery_present, original_input,
        recovery_input, restored);
    require(restored_word == original[missing],
        "maximum high R=1 recovery mismatch");
    compared_bytes += 2;
    leo2_codec_destroy(codec);
}

void test_high_general_limit(leo2_context* context, uint64_t& compared_bytes)
{
    const uint32_t k = 65534;
    const uint32_t r = 2;
    leo2_codec* codec = make_codec(
        context, k, r, LEO2_PROFILE_LEGACY_HIGH_V1);
    require(leo2_codec_padded_side(codec) == 2,
        "maximum high general padded side mismatch");

    std::vector<Word> original(k);
    std::vector<const void*> original_input(k);
    for (uint32_t i = 0; i < k; ++i)
    {
        original[i][0] = static_cast<unsigned char>(i * 11u + 0x51u);
        original[i][1] = static_cast<unsigned char>(i * 23u + (i >> 7));
        original_input[i] = original[i].data();
    }
    std::array<Word, 2> parity = {{{{ 0, 0 }}, {{ 0, 0 }}}};
    std::vector<void*> parity_output(r);
    for (uint32_t i = 0; i < r; ++i)
        parity_output[i] = parity[i].data();
    AlignedBuffer scratch = encode_scratch(codec);
    require_result(leo2_encode(codec, 2, &original_input[0],
        &parity_output[0], scratch.data(), scratch.size()),
        "maximum high general encode");

    const uint32_t first_missing = 0;
    const uint32_t second_missing = k - 1;
    std::vector<uint8_t> original_present(k, 1);
    original_present[first_missing] = 0;
    original_present[second_missing] = 0;
    std::vector<uint8_t> recovery_present(r, 1);
    original_input[first_missing] = NULL;
    original_input[second_missing] = NULL;
    std::vector<const void*> recovery_input(r);
    for (uint32_t i = 0; i < r; ++i)
        recovery_input[i] = parity[i].data();
    std::array<Word, 2> restored_words = {{{{ 0, 0 }}, {{ 0, 0 }}}};
    std::vector<void*> restored(k, NULL);
    restored[first_missing] = restored_words[0].data();
    restored[second_missing] = restored_words[1].data();
    recover_missing(codec, original_present, recovery_present, original_input,
        recovery_input, restored);
    require(restored_words[0] == original[first_missing] &&
        restored_words[1] == original[second_missing],
        "maximum high general recovery mismatch");
    compared_bytes += 4;
    leo2_codec_destroy(codec);
}

void test_low_limit(leo2_context* context, uint64_t& compared_bytes)
{
    const uint32_t k = 2;
    const uint32_t r = 65534;
    leo2_codec* codec = make_codec(context, k, r, LEO2_PROFILE_LOW_V1);
    require(leo2_codec_padded_side(codec) == 2,
        "maximum low padded side mismatch");

    std::array<Word, 2> original = {{{{ 0x5a, 0xc3 }}, {{ 0x17, 0xe1 }}}};
    std::vector<const void*> original_input(k);
    for (uint32_t i = 0; i < k; ++i)
        original_input[i] = original[i].data();

    std::array<Word, 2> parity = {{{{ 0, 0 }}, {{ 0, 0 }}}};
    std::vector<void*> parity_output(r, NULL);
    parity_output[r - 2] = parity[0].data();
    parity_output[r - 1] = parity[1].data();
    AlignedBuffer scratch = encode_scratch(codec);
    require_result(leo2_encode(codec, 2, &original_input[0],
        &parity_output[0], scratch.data(), scratch.size()),
        "maximum low encode");

    std::vector<uint8_t> original_present(k, 0);
    std::vector<uint8_t> recovery_present(r, 0);
    recovery_present[r - 2] = 1;
    recovery_present[r - 1] = 1;
    original_input[0] = NULL;
    original_input[1] = NULL;
    std::vector<const void*> recovery_input(r, NULL);
    recovery_input[r - 2] = parity[0].data();
    recovery_input[r - 1] = parity[1].data();
    std::array<Word, 2> restored_words = {{{{ 0, 0 }}, {{ 0, 0 }}}};
    std::vector<void*> restored(k, NULL);
    restored[0] = restored_words[0].data();
    restored[1] = restored_words[1].data();
    recover_missing(codec, original_present, recovery_present, original_input,
        recovery_input, restored);
    require(restored_words == original, "maximum low recovery mismatch");
    compared_bytes += 4;
    leo2_codec_destroy(codec);
}

void test_low_k1_literal_limit(
    leo2_context* context,
    leo2_field field,
    uint32_t r,
    uint32_t expected_parent,
    uint64_t& compared_bytes)
{
    leo2_codec* codec = NULL;
    require_result(leo2_codec_create(context, 1, r, LEO2_PROFILE_LOW_V1,
        field, NULL, &codec), "literal maximum low codec create");
    require(codec != NULL, "literal maximum low codec create returned null");
    require(leo2_codec_parent_count(codec) == expected_parent,
        "literal maximum low parent mismatch");
    require(leo2_codec_padded_side(codec) == 1,
        "literal maximum low padded side mismatch");

    const Word original = {{ 0x6d, 0xb2 }};
    std::vector<const void*> original_input(1, original.data());
    std::array<Word, 2> parity = {{{{ 0, 0 }}, {{ 0, 0 }}}};
    std::vector<void*> parity_output(r, NULL);
    parity_output[0] = parity[0].data();
    parity_output[r - 1] = parity[1].data();
    AlignedBuffer scratch = encode_scratch(codec);
    require_result(leo2_encode(codec, 2, &original_input[0],
        &parity_output[0], scratch.data(), scratch.size()),
        "literal maximum low encode");
    require(parity[0] == original,
        "literal maximum low first parity is not the sole original");
    require(parity[1] == original,
        "literal maximum low last parity is not the sole original");
    compared_bytes += 4;

    std::vector<uint8_t> original_present(1, 0);
    std::vector<uint8_t> recovery_present(r, 0);
    recovery_present[r - 1] = 1;
    original_input[0] = NULL;
    std::vector<const void*> recovery_input(r, NULL);
    recovery_input[r - 1] = parity[1].data();
    Word restored_word = {{ 0, 0 }};
    std::vector<void*> restored(1, restored_word.data());
    recover_missing(codec, original_present, recovery_present, original_input,
        recovery_input, restored);
    require(restored_word == original,
        "literal maximum low single-parity recovery mismatch");
    compared_bytes += 2;
    leo2_codec_destroy(codec);
}

void test_limit_rejection(leo2_context* context)
{
    leo2_codec* codec = NULL;
    require(leo2_codec_create(context, 1, 256,
        LEO2_PROFILE_LOW_V1, LEO2_FIELD_GF8, NULL, &codec) ==
        LEO2_INVALID_COUNTS, "GF8 low parent overflow was accepted");
    require(codec == NULL, "GF8 low overflow returned a codec");
    require(leo2_codec_create(context, 65535, 2,
        LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF16, NULL, &codec) ==
        LEO2_INVALID_COUNTS, "GF16 high parent overflow was accepted");
    require(codec == NULL, "GF16 high overflow returned a codec");
    require(leo2_codec_create(context, 1, 65536,
        LEO2_PROFILE_LOW_V1, LEO2_FIELD_GF16, NULL, &codec) ==
        LEO2_INVALID_COUNTS, "GF16 low parent overflow was accepted");
    require(codec == NULL, "GF16 low overflow returned a codec");
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
        require_result(leo2_context_create(&options, &context),
            "maximum context create");

        uint64_t compared_bytes = 0;
        test_high_r1_limit(context, compared_bytes);
        test_high_general_limit(context, compared_bytes);
        test_low_limit(context, compared_bytes);
        test_low_k1_literal_limit(context, LEO2_FIELD_GF8, 255, 256,
            compared_bytes);
        test_low_k1_literal_limit(context, LEO2_FIELD_GF16, 65535, 65536,
            compared_bytes);
        test_limit_rejection(context);

        leo2_context_destroy(context);
        std::cout << "leopard2 maximum-count test passed: full_parents=5"
                  << " literal_low_k1_profiles=2"
                  << " compared_bytes=" << compared_bytes << std::endl;
        return 0;
    }
    catch (const std::exception& error)
    {
        std::cerr << "leopard2 maximum-count test failed: "
                  << error.what() << std::endl;
        return 1;
    }
}
