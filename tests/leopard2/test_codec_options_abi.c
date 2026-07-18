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
    struct short_codec_v1 {
        size_t struct_size;
        uint32_t flags;
        uint32_t reserved;
        leo2_codec *output;
    };
    union short_codec_options_storage {
        leo2_codec_options current;
        struct short_codec_v1 v1;
    } short_options;
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

    {
        leo2_context_options aliased_options = context_options;
        leo2_context_options snapshot = aliased_options;
        if (!require_result(leo2_context_create(&aliased_options,
                (leo2_context **)(void *)&aliased_options), LEO2_OVERLAP,
                "C ABI context options/output overlap") ||
            memcmp(&aliased_options, &snapshot, sizeof(snapshot)) != 0) {
            leo2_context_destroy(context);
            return 1;
        }
    }
    {
        leo2_context *untouched = (leo2_context *)(uintptr_t)1;
        leo2_context_options impossible_options = context_options;
        impossible_options.struct_size = SIZE_MAX;
        if (!require_result(leo2_context_create(&impossible_options,
                &untouched), LEO2_INVALID_ARGUMENT,
                "C ABI unrepresentable context options span") ||
            untouched != (leo2_context *)(uintptr_t)1) {
            leo2_context_destroy(context);
            return 1;
        }
        if (!require_result(leo2_context_create(
                (const leo2_context_options *)(uintptr_t)(UINTPTR_MAX -
                    sizeof(size_t) + 1), &untouched),
                LEO2_INVALID_ARGUMENT,
                "C ABI unrepresentable context options prefix") ||
            untouched != (leo2_context *)(uintptr_t)1) {
            leo2_context_destroy(context);
            return 1;
        }
        if (!require_result(leo2_context_create(NULL,
                (leo2_context **)(uintptr_t)(UINTPTR_MAX -
                    sizeof(leo2_context *) + 1)), LEO2_INVALID_ARGUMENT,
                "C ABI unrepresentable context output span")) {
            leo2_context_destroy(context);
            return 1;
        }
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

    {
        leo2_codec_options aliased_options = options;
        leo2_codec_options snapshot = aliased_options;
        if (!require_result(leo2_codec_create(context, 5, 3,
                LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF16,
                &aliased_options,
                (leo2_codec **)(void *)&aliased_options), LEO2_OVERLAP,
                "C ABI codec options/output overlap") ||
            memcmp(&aliased_options, &snapshot, sizeof(snapshot)) != 0) {
            leo2_context_destroy(context);
            return 1;
        }
    }

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

    {
        struct extended_codec_options snapshot = extended_options;
        if (!require_result(leo2_codec_create(context, 5, 3,
                LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8,
                &extended_options.base,
                (leo2_codec **)(void *)&extended_options.future_words[0]),
                LEO2_OVERLAP,
                "C ABI extension/output overlap") ||
            memcmp(&extended_options, &snapshot, sizeof(snapshot)) != 0) {
            leo2_context_destroy(context);
            return 1;
        }
    }

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

    /* The output immediately after an exact v1 prefix is not part of the
       declared immutable input and remains a valid forward-compatible use. */
    memset(&short_options, 0, sizeof(short_options));
    if (offsetof(struct short_codec_v1, output) != version1_size) {
        leo2_context_destroy(context);
        return 1;
    }
    short_options.v1.struct_size = version1_size;
    if (!require_result(leo2_codec_create(context, 5, 3,
            LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF16,
            &short_options.current, &short_options.v1.output),
            LEO2_SUCCESS, "C ABI adjacent v1-prefix output") ||
        short_options.v1.output == NULL ||
        leo2_codec_shard_layout(short_options.v1.output) !=
            LEO2_SHARD_LAYOUT_NATIVE_V1) {
        leo2_codec_destroy(short_options.v1.output);
        leo2_context_destroy(context);
        return 1;
    }
    leo2_codec_destroy(short_options.v1.output);
    short_options.v1.output = NULL;

    memset(&options, 0, sizeof(options));
    options.struct_size = SIZE_MAX;
    codec = (leo2_codec *)(uintptr_t)1;
    if (!require_result(leo2_codec_create(context, 5, 3,
            LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8,
            &options, &codec), LEO2_INVALID_ARGUMENT,
            "C ABI unrepresentable codec options span") ||
        codec != (leo2_codec *)(uintptr_t)1) {
        leo2_context_destroy(context);
        return 1;
    }
    if (!require_result(leo2_codec_create(context, 5, 3,
            LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8,
            (const leo2_codec_options *)(uintptr_t)(UINTPTR_MAX -
                sizeof(size_t) + 1), &codec), LEO2_INVALID_ARGUMENT,
            "C ABI unrepresentable codec options prefix") ||
        codec != (leo2_codec *)(uintptr_t)1 ||
        !require_result(leo2_codec_create(context, 5, 3,
            LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8, NULL,
            (leo2_codec **)(uintptr_t)(UINTPTR_MAX -
                sizeof(leo2_codec *) + 1)), LEO2_INVALID_ARGUMENT,
            "C ABI unrepresentable codec output span")) {
        leo2_context_destroy(context);
        return 1;
    }

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

    {
        leo2_codec *plan_codec = NULL;
        uint8_t original_present[9];
        uint8_t recovery_present[7];
        union original_presence_alias {
            void *alignment;
            uint8_t bytes[9];
        } original_alias;
        union recovery_presence_alias {
            void *alignment;
            uint8_t bytes[sizeof(void *) > 7 ? sizeof(void *) : 7];
        } recovery_alias;
        uint8_t original_snapshot[sizeof(original_alias.bytes)];
        uint8_t recovery_snapshot[sizeof(recovery_alias.bytes)];
        leo2_decode_plan *plan = NULL;

        memset(original_present, 1, sizeof(original_present));
        memset(recovery_present, 1, sizeof(recovery_present));
        if (!require_result(leo2_codec_create(context, 9, 7,
                LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8, NULL,
                &plan_codec), LEO2_SUCCESS,
                "C ABI presence-alias codec create")) {
            leo2_context_destroy(context);
            return 1;
        }

        memset(original_alias.bytes, 1, sizeof(original_alias.bytes));
        memcpy(original_snapshot, original_alias.bytes,
            sizeof(original_snapshot));
        if (!require_result(leo2_decode_plan_create(plan_codec,
                original_alias.bytes, recovery_present,
                (leo2_decode_plan **)(void *)original_alias.bytes),
                LEO2_OVERLAP, "C ABI original-presence/output overlap") ||
            memcmp(original_alias.bytes, original_snapshot,
                sizeof(original_snapshot)) != 0) {
            leo2_codec_destroy(plan_codec);
            leo2_context_destroy(context);
            return 1;
        }

        memset(original_alias.bytes, 1, sizeof(original_alias.bytes));
        memcpy(original_snapshot, original_alias.bytes,
            sizeof(original_snapshot));
        if (!require_result(leo2_decode_plan_create(plan_codec,
                original_alias.bytes, NULL,
                (leo2_decode_plan **)(void *)original_alias.bytes),
                LEO2_OVERLAP,
                "C ABI original-presence/output overlap with null recovery") ||
            memcmp(original_alias.bytes, original_snapshot,
                sizeof(original_snapshot)) != 0) {
            leo2_codec_destroy(plan_codec);
            leo2_context_destroy(context);
            return 1;
        }

        memset(recovery_alias.bytes, 0xa5, sizeof(recovery_alias.bytes));
        memset(recovery_alias.bytes, 1, sizeof(recovery_present));
        memcpy(recovery_snapshot, recovery_alias.bytes,
            sizeof(recovery_snapshot));
        if (!require_result(leo2_decode_plan_create(plan_codec,
                original_present, recovery_alias.bytes,
                (leo2_decode_plan **)(void *)recovery_alias.bytes),
                LEO2_OVERLAP, "C ABI recovery-presence/output overlap") ||
            memcmp(recovery_alias.bytes, recovery_snapshot,
                sizeof(recovery_snapshot)) != 0) {
            leo2_codec_destroy(plan_codec);
            leo2_context_destroy(context);
            return 1;
        }

        memset(recovery_alias.bytes, 0xa5, sizeof(recovery_alias.bytes));
        memset(recovery_alias.bytes, 1, sizeof(recovery_present));
        memcpy(recovery_snapshot, recovery_alias.bytes,
            sizeof(recovery_snapshot));
        if (!require_result(leo2_decode_plan_create(plan_codec,
                NULL, recovery_alias.bytes,
                (leo2_decode_plan **)(void *)recovery_alias.bytes),
                LEO2_OVERLAP,
                "C ABI recovery-presence/output overlap with null original") ||
            memcmp(recovery_alias.bytes, recovery_snapshot,
                sizeof(recovery_snapshot)) != 0) {
            leo2_codec_destroy(plan_codec);
            leo2_context_destroy(context);
            return 1;
        }

        plan = (leo2_decode_plan *)(uintptr_t)1;
        if (!require_result(leo2_decode_plan_create(plan_codec,
                (const uint8_t *)(uintptr_t)(UINTPTR_MAX - 4),
                recovery_present, &plan), LEO2_INVALID_ARGUMENT,
                "C ABI unrepresentable original-presence span") ||
            plan != (leo2_decode_plan *)(uintptr_t)1) {
            leo2_codec_destroy(plan_codec);
            leo2_context_destroy(context);
            return 1;
        }
        plan = (leo2_decode_plan *)(uintptr_t)1;
        if (!require_result(leo2_decode_plan_create(plan_codec,
                original_present,
                (const uint8_t *)(uintptr_t)(UINTPTR_MAX - 3),
                &plan), LEO2_INVALID_ARGUMENT,
                "C ABI unrepresentable recovery-presence span") ||
            plan != (leo2_decode_plan *)(uintptr_t)1 ||
            !require_result(leo2_decode_plan_create(plan_codec,
                original_present, recovery_present,
                (leo2_decode_plan **)(uintptr_t)(UINTPTR_MAX -
                    sizeof(leo2_decode_plan *) + 1)),
                LEO2_INVALID_ARGUMENT,
                "C ABI unrepresentable plan output span")) {
            leo2_codec_destroy(plan_codec);
            leo2_context_destroy(context);
            return 1;
        }

        original_present[0] = 2;
        plan = (leo2_decode_plan *)(uintptr_t)1;
        if (!require_result(leo2_decode_plan_create(plan_codec,
                original_present, recovery_present, &plan),
                LEO2_INVALID_ARGUMENT,
                "C ABI ordinary invalid presence rejection") ||
            plan != NULL) {
            leo2_codec_destroy(plan_codec);
            leo2_context_destroy(context);
            return 1;
        }
        leo2_codec_destroy(plan_codec);
    }

#if (defined(__unix__) || defined(__APPLE__)) && defined(MAP_ANONYMOUS)
    {
        const long page_size_value = sysconf(_SC_PAGESIZE);
        leo2_codec *guard_codec = NULL;
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

        if (!require_result(leo2_codec_create(context, 16, 8,
                LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8, NULL,
                &guard_codec), LEO2_SUCCESS,
                "C ABI guard-presence codec create")) {
            munmap(mapping, page_size * 2);
            leo2_context_destroy(context);
            return 1;
        }

        /* Only struct_size is mapped.  Reading flags before checking this
           deliberately short prefix crosses into the protected page. */
        leo2_codec_options *guarded_options = (leo2_codec_options *)(
            mapping + page_size - sizeof(size_t));
        guarded_options->struct_size = sizeof(size_t);
        {
            leo2_context *guarded_context =
                (leo2_context *)(uintptr_t)1;
            if (!require_result(leo2_context_create(
                    (const leo2_context_options *)guarded_options,
                    &guarded_context), LEO2_INVALID_ARGUMENT,
                    "C ABI guard-page context prefix rejection") ||
                guarded_context != NULL) {
                leo2_codec_destroy(guard_codec);
                munmap(mapping, page_size * 2);
                leo2_context_destroy(context);
                return 1;
            }
        }
        codec = (leo2_codec *)(uintptr_t)1;
        if (!require_result(leo2_codec_create(context, 5, 3,
                LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8,
                guarded_options, &codec), LEO2_INVALID_ARGUMENT,
                "C ABI guard-page prefix rejection") || codec != NULL) {
            leo2_codec_destroy(guard_codec);
            munmap(mapping, page_size * 2);
            leo2_codec_destroy(codec);
            leo2_context_destroy(context);
            return 1;
        }

        {
            uint8_t *guarded_original = mapping + page_size - 16;
            uint8_t *guarded_recovery = guarded_original - 8;
            uint8_t original_snapshot[16];
            leo2_decode_plan *guarded_plan = NULL;
            memset(guarded_original, 1, 16);
            memset(guarded_recovery, 1, 8);
            if (!require_result(leo2_decode_plan_create(guard_codec,
                    guarded_original, guarded_recovery, &guarded_plan),
                    LEO2_SUCCESS,
                    "C ABI exact guard-page presence spans")) {
                leo2_codec_destroy(guard_codec);
                munmap(mapping, page_size * 2);
                leo2_context_destroy(context);
                return 1;
            }
            leo2_decode_plan_destroy(guarded_plan);
            memcpy(original_snapshot, guarded_original,
                sizeof(original_snapshot));
            if (!require_result(leo2_decode_plan_create(guard_codec,
                    guarded_original, guarded_recovery,
                    (leo2_decode_plan **)(void *)guarded_original),
                    LEO2_OVERLAP,
                    "C ABI guard-page presence/output overlap") ||
                memcmp(guarded_original, original_snapshot,
                    sizeof(original_snapshot)) != 0) {
                leo2_codec_destroy(guard_codec);
                munmap(mapping, page_size * 2);
                leo2_context_destroy(context);
                return 1;
            }
        }
        leo2_codec_destroy(guard_codec);
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
