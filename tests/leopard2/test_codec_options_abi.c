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

#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

#if defined(__unix__) || defined(__APPLE__)
#include <sys/mman.h>
#include <unistd.h>
#endif

typedef char leo2_layout_storage_must_be_u32[
    sizeof(((leo2_codec_options *)0)->shard_layout) == sizeof(uint32_t) ? 1 : -1];
typedef char leo2_backend_storage_must_be_u32[
    sizeof(((leo2_context_options *)0)->backend) == sizeof(uint32_t) ? 1 : -1];

static int require_result(leo2_result actual, leo2_result expected,
    const char *operation)
{
    if (actual == expected) {
        return 1;
    }
    fprintf(stderr, "%s: got %s (%d), expected %s (%d)\n", operation,
        leo2_result_string(actual), (int)actual,
        leo2_result_string(expected), (int)expected);
    return 0;
}

int main(void)
{
    leo2_context *context = NULL;
    leo2_codec *codec = NULL;
    leo2_context_options context_options;
    leo2_codec_options options;
    struct extended_codec_options {
        leo2_codec_options base;
        uint64_t future_words[2];
    } extended_options;
    const size_t version1_size = offsetof(leo2_codec_options, shard_layout);
    const size_t layout_field_end = version1_size + sizeof(options.shard_layout);

    memset(&context_options, 0, sizeof(context_options));
    context_options.struct_size = sizeof(context_options);
    context_options.thread_count = 1;
    if (LEO2_API_VERSION < 2u ||
        !require_result(leo2_context_create(&context_options, &context), LEO2_SUCCESS,
            "C ABI context create")) {
        return 1;
    }

    /* Leave trailing struct padding nonzero.  A short-enum C caller must still
       initialize the complete fixed-width layout field seen by the library. */
    memset(&options, 0xa5, sizeof(options));
    options.struct_size = sizeof(options);
    options.flags = 0;
    options.reserved = 0;
    options.shard_layout = LEO2_SHARD_LAYOUT_GF16_PADDED_ODD_V1;
    if (!require_result(leo2_codec_create(context, 5, 3,
            LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF16, &options, &codec),
            LEO2_SUCCESS, "C ABI padded codec create") ||
        leo2_codec_shard_layout(codec) !=
            LEO2_SHARD_LAYOUT_GF16_PADDED_ODD_V1) {
        leo2_codec_destroy(codec);
        leo2_context_destroy(context);
        return 1;
    }
    leo2_codec_destroy(codec);
    codec = NULL;

    /* A future caller may append fields.  The known current prefix remains
       valid and unknown trailing storage is ignored. */
    memset(&extended_options, 0xa5, sizeof(extended_options));
    extended_options.base.struct_size = sizeof(extended_options);
    extended_options.base.flags = 0;
    extended_options.base.reserved = 0;
    extended_options.base.shard_layout = LEO2_SHARD_LAYOUT_NATIVE_V1;
    if (!require_result(leo2_codec_create(context, 5, 3,
            LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8,
            &extended_options.base, &codec), LEO2_SUCCESS,
            "C ABI oversized-extension codec create") ||
        leo2_codec_shard_layout(codec) != LEO2_SHARD_LAYOUT_NATIVE_V1) {
        leo2_codec_destroy(codec);
        leo2_context_destroy(context);
        return 1;
    }
    leo2_codec_destroy(codec);
    codec = NULL;

    memset(&options, 0, sizeof(options));
    options.struct_size = sizeof(options);
    options.reserved = 1;
    codec = (leo2_codec *)(uintptr_t)1;
    if (!require_result(leo2_codec_create(context, 5, 3,
            LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8, &options, &codec),
            LEO2_INVALID_ARGUMENT, "C ABI codec reserved rejection") ||
        codec != NULL) {
        leo2_context_destroy(context);
        return 1;
    }

    /* The exact v1 prefix ignores bytes belonging to the appended field. */
    memset(&options, 0xa5, sizeof(options));
    options.struct_size = version1_size;
    options.flags = 0;
    options.reserved = 0;
    if (!require_result(leo2_codec_create(context, 5, 3,
            LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF16, &options, &codec),
            LEO2_SUCCESS, "C ABI v1-prefix codec create") ||
        leo2_codec_shard_layout(codec) != LEO2_SHARD_LAYOUT_NATIVE_V1) {
        leo2_codec_destroy(codec);
        leo2_context_destroy(context);
        return 1;
    }
    leo2_codec_destroy(codec);
    codec = NULL;

    memset(&options, 0, sizeof(options));
    options.struct_size = layout_field_end - 1;
    if (!require_result(leo2_codec_create(context, 5, 3,
            LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF16, &options, &codec),
            LEO2_INVALID_ARGUMENT, "C ABI partial-field rejection") ||
        codec != NULL) {
        leo2_codec_destroy(codec);
        leo2_context_destroy(context);
        return 1;
    }

    /* An undersized portable object is rejected before any later field is
       interpreted, regardless of the bytes present in that storage. */
    memset(&options, 0xa5, sizeof(options));
    options.struct_size = version1_size - 1;
    codec = (leo2_codec *)(uintptr_t)1;
    if (!require_result(leo2_codec_create(context, 5, 3,
            LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8, &options, &codec),
            LEO2_INVALID_ARGUMENT, "C ABI undersized-prefix rejection") ||
        codec != NULL) {
        leo2_codec_destroy(codec);
        leo2_context_destroy(context);
        return 1;
    }

#if (defined(__unix__) || defined(__APPLE__)) && defined(MAP_ANONYMOUS)
    {
        const long page_size_value = sysconf(_SC_PAGESIZE);
        if (page_size_value <= 0) {
            leo2_context_destroy(context);
            return 1;
        }
        const size_t page_size = (size_t)page_size_value;
        unsigned char *mapping = (unsigned char *)mmap(NULL, page_size * 2,
            PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
        if (mapping == MAP_FAILED ||
            mprotect(mapping + page_size, page_size, PROT_NONE) != 0) {
            if (mapping != MAP_FAILED)
                munmap(mapping, page_size * 2);
            leo2_context_destroy(context);
            return 1;
        }

        /* Only struct_size is mapped.  Reading flags before checking this
           deliberately short prefix crosses into the protected page. */
        leo2_codec_options *guarded_options = (leo2_codec_options *)(
            mapping + page_size - sizeof(size_t));
        guarded_options->struct_size = sizeof(size_t);
        codec = (leo2_codec *)(uintptr_t)1;
        if (!require_result(leo2_codec_create(context, 5, 3,
                LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8,
                guarded_options, &codec), LEO2_INVALID_ARGUMENT,
                "C ABI guard-page prefix rejection") || codec != NULL) {
            munmap(mapping, page_size * 2);
            leo2_codec_destroy(codec);
            leo2_context_destroy(context);
            return 1;
        }
        if (munmap(mapping, page_size * 2) != 0) {
            leo2_context_destroy(context);
            return 1;
        }
    }
#endif

    leo2_context_destroy(context);
    context = (leo2_context *)(uintptr_t)1;
    memset(&context_options, 0, sizeof(context_options));
    context_options.struct_size = sizeof(context_options);
    context_options.reserved = 1;
    if (!require_result(leo2_context_create(&context_options, &context),
            LEO2_INVALID_ARGUMENT, "C ABI context reserved rejection") ||
        context != NULL) {
        return 1;
    }

    printf("leopard2 C options ABI passed: backend_bytes=%lu layout_bytes=%lu v1_bytes=%lu\n",
        (unsigned long)sizeof(context_options.backend),
        (unsigned long)sizeof(options.shard_layout),
        (unsigned long)version1_size);
    return 0;
}
