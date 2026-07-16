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

#ifndef CAT_LEOPARD2_RS_H
#define CAT_LEOPARD2_RS_H

#include <stddef.h>
#include <stdint.h>

#if defined(LEO_BUILDING)
# if defined(LEO_DLL)
#  define LEO2_EXPORT __declspec(dllexport)
# else
#  define LEO2_EXPORT
# endif
#else
# if defined(LEO_DLL)
#  define LEO2_EXPORT __declspec(dllimport)
# else
#  define LEO2_EXPORT extern
# endif
#endif

#ifdef __cplusplus
extern "C" {
#endif

#define LEO2_API_VERSION 1u

typedef struct leo2_context leo2_context;
typedef struct leo2_codec leo2_codec;
typedef struct leo2_decode_plan leo2_decode_plan;

typedef enum leo2_result {
    LEO2_SUCCESS = 0,
    LEO2_NEED_MORE_DATA = -1,
    LEO2_INVALID_ARGUMENT = -2,
    LEO2_INVALID_COUNTS = -3,
    LEO2_UNSUPPORTED = -4,
    LEO2_SCRATCH_TOO_SMALL = -5,
    LEO2_BAD_ALIGNMENT = -6,
    LEO2_OVERLAP = -7,
    LEO2_OUT_OF_MEMORY = -8,
    LEO2_INTERNAL_ERROR = -9
} leo2_result;

typedef enum leo2_profile {
    LEO2_PROFILE_AUTO = 0,
    LEO2_PROFILE_LEGACY_HIGH_V1 = 1,
    LEO2_PROFILE_LOW_V1 = 2,
    LEO2_PROFILE_EXACT_EXPERIMENTAL_V1 = 3
} leo2_profile;

typedef enum leo2_field {
    LEO2_FIELD_AUTO = 0,
    LEO2_FIELD_GF8 = 1,
    LEO2_FIELD_GF16 = 2
} leo2_field;

typedef enum leo2_backend {
    LEO2_BACKEND_AUTO = 0,
    LEO2_BACKEND_SCALAR = 1,
    LEO2_BACKEND_SSSE3 = 2,
    LEO2_BACKEND_AVX2 = 3,
    LEO2_BACKEND_NEON = 4
} leo2_backend;

/*
    A future caller-executor adapter can be added through a larger struct_size.
    thread_count selects the persistent context pool used by batch entry points;
    the calling thread participates and no worker is created per batch call.
*/
typedef struct leo2_context_options {
    size_t struct_size;
    leo2_backend backend;
    uint32_t thread_count; /* zero means library default */
    uint32_t reserved;
} leo2_context_options;

typedef struct leo2_codec_options {
    size_t struct_size;
    uint32_t flags;
    uint32_t reserved;
} leo2_codec_options;

/* Test/diagnostic flag retaining the full O(N log N) decoder fallback. */
#define LEO2_CODEC_FORCE_GENERIC_DECODE 0x00000001u

typedef struct leo2_encode_batch_item {
    uint64_t shard_bytes;
    const void* const* original;
    void* const* recovery;
    void* scratch;
    size_t scratch_bytes;
} leo2_encode_batch_item;

typedef struct leo2_decode_batch_item {
    uint64_t shard_bytes;
    const void* const* original;
    const void* const* recovery;
    void* const* restored_original;
    void* scratch;
    size_t scratch_bytes;
} leo2_decode_batch_item;

LEO2_EXPORT const char* leo2_result_string(leo2_result result);

LEO2_EXPORT leo2_result leo2_context_create(
    const leo2_context_options* options,
    leo2_context** context_out);
LEO2_EXPORT void leo2_context_destroy(leo2_context* context);
LEO2_EXPORT leo2_backend leo2_context_backend(const leo2_context* context);
LEO2_EXPORT uint32_t leo2_context_thread_count(const leo2_context* context);

LEO2_EXPORT leo2_result leo2_codec_create(
    leo2_context* context,
    uint32_t original_count,
    uint32_t recovery_count,
    leo2_profile profile,
    leo2_field field,
    const leo2_codec_options* options,
    leo2_codec** codec_out);
LEO2_EXPORT void leo2_codec_destroy(leo2_codec* codec);

LEO2_EXPORT uint32_t leo2_codec_original_count(const leo2_codec* codec);
LEO2_EXPORT uint32_t leo2_codec_recovery_count(const leo2_codec* codec);
LEO2_EXPORT uint32_t leo2_codec_parent_count(const leo2_codec* codec);
LEO2_EXPORT uint32_t leo2_codec_padded_side(const leo2_codec* codec);
LEO2_EXPORT leo2_profile leo2_codec_profile(const leo2_codec* codec);
LEO2_EXPORT leo2_field leo2_codec_field(const leo2_codec* codec);

/*
    Scratch must start at leo2_scratch_alignment() alignment.  Scratch and all
    non-null shard ranges must be mutually disjoint, except that input shards may
    alias other input shards.  Recovery output entries may be null to request a
    parity subset.  Encode/decode execution performs no allocation.
*/
LEO2_EXPORT size_t leo2_scratch_alignment(void);
LEO2_EXPORT leo2_result leo2_encode_scratch_size(
    const leo2_codec* codec,
    uint64_t shard_bytes,
    size_t* scratch_bytes_out);
LEO2_EXPORT leo2_result leo2_encode(
    const leo2_codec* codec,
    uint64_t shard_bytes,
    const void* const* original,
    void* const* recovery,
    void* scratch,
    size_t scratch_bytes);
LEO2_EXPORT leo2_result leo2_encode_batch(
    const leo2_codec* codec,
    const leo2_encode_batch_item* items,
    size_t item_count);

/*
    Presence arrays contain one for a received shard and zero for an erasure.
    They are copied into the immutable plan.  The codec and its context must
    outlive the plan.  Execution restores missing originals only; parity rebuild
    is explicit encoding after the caller combines surviving/restored originals.
*/
LEO2_EXPORT leo2_result leo2_decode_plan_create(
    const leo2_codec* codec,
    const uint8_t* original_present,
    const uint8_t* recovery_present,
    leo2_decode_plan** plan_out);
LEO2_EXPORT void leo2_decode_plan_destroy(leo2_decode_plan* plan);
LEO2_EXPORT uint32_t leo2_decode_plan_missing_original_count(
    const leo2_decode_plan* plan);
LEO2_EXPORT leo2_result leo2_decode_plan_scratch_size(
    const leo2_decode_plan* plan,
    uint64_t shard_bytes,
    size_t* scratch_bytes_out);
LEO2_EXPORT leo2_result leo2_decode_plan_execute(
    const leo2_decode_plan* plan,
    uint64_t shard_bytes,
    const void* const* original,
    const void* const* recovery,
    void* const* restored_original,
    void* scratch,
    size_t scratch_bytes);
LEO2_EXPORT leo2_result leo2_decode_plan_execute_batch(
    const leo2_decode_plan* plan,
    const leo2_decode_batch_item* items,
    size_t item_count);

/* One-shot convenience wrapper.  It allocates and destroys plan setup state. */
LEO2_EXPORT leo2_result leo2_decode_scratch_size(
    const leo2_codec* codec,
    uint64_t shard_bytes,
    size_t* scratch_bytes_out);
LEO2_EXPORT leo2_result leo2_decode(
    const leo2_codec* codec,
    uint64_t shard_bytes,
    const uint8_t* original_present,
    const uint8_t* recovery_present,
    const void* const* original,
    const void* const* recovery,
    void* const* restored_original,
    void* scratch,
    size_t scratch_bytes);

#ifdef __cplusplus
}
#endif

#endif /* CAT_LEOPARD2_RS_H */
