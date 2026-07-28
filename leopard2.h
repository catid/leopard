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

#define LEO2_API_VERSION 4u

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

/*
    EXACT_EXPERIMENTAL_V1 reserves a distinct code identity for research.  It is
    not implemented by the stable codec constructors and currently returns
    LEO2_UNSUPPORTED; it never aliases either production V1 profile.
*/

typedef enum leo2_field {
    LEO2_FIELD_AUTO = 0,
    LEO2_FIELD_GF8 = 1,
    LEO2_FIELD_GF16 = 2
} leo2_field;

/* Bits returned by leo2_context_field_mask(). */
#define LEO2_FIELD_MASK_GF8  0x00000001u
#define LEO2_FIELD_MASK_GF16 0x00000002u

/*
    Reports mathematical field implementations usable by this initialized
    context.  The value is independent of the selected CPU backend.
    AUTO field selection remains canonical and build-independent: a code whose
    canonical field is omitted is unsupported rather than silently changing
    its wire representation.
*/
LEO2_EXPORT uint32_t leo2_context_field_mask(const leo2_context* context);

/*
    Shard layout is part of persistent code identity.  Native V1 retains the
    existing GF8 and even-byte GF16 representation.  Padded-odd V1 is an
    explicit GF16-only framing profile: an odd B-byte application payload is
    stored as a physical W=B+1 byte systematic shard with wire[B] equal to
    zero.  Every parity shard also stores all W bytes; its last byte is not
    padding and may be nonzero.
*/
typedef enum leo2_shard_layout {
    LEO2_SHARD_LAYOUT_NATIVE_V1 = 0,
    LEO2_SHARD_LAYOUT_GF16_PADDED_ODD_V1 = 1
} leo2_shard_layout;

typedef enum leo2_backend {
    LEO2_BACKEND_AUTO = 0,
    LEO2_BACKEND_SCALAR = 1,
    LEO2_BACKEND_SSSE3 = 2,
    LEO2_BACKEND_AVX2 = 3,
    LEO2_BACKEND_NEON = 4,
    /* Explicit AVX-512F/BW/VL acceleration. */
    LEO2_BACKEND_AVX512 = 5,
    /* Explicit AVX2 + GFNI acceleration (VEX-encoded VGF2P8AFFINEQB). */
    LEO2_BACKEND_GFNI = 6
} leo2_backend;

/*
    Option structs must be zero-initialized before their documented fields are
    assigned.  reserved must be zero.  struct_size may exceed the current
    definition; the known prefix is consumed and an unknown trailing extension
    is ignored.  A future caller-executor adapter can therefore be added without
    changing the entry point.

    thread_count selects the persistent context pool used by batch entry points;
    the calling thread participates.  Workers are started lazily by the first
    batch containing more than one item, grown only when a larger batch needs
    more parallelism, then reused without recreating existing workers.
    A zero value uses the process's allowed CPU count where the platform exposes
    it, capped at 128; explicit nonzero values retain their requested budget.
    An explicit backend selects one immutable execution table for this context.
    AUTO uses the production-default startup-qualified runtime backend as its
    reported baseline.  On production
    x86 builds, SCALAR, SSSE3, AVX2, AVX512, and GFNI may explicitly
    select any compiled, host-qualified table.  AUTO deliberately reports AVX2
    as its baseline because AVX512 whole-codec evidence supports only bounded
    encode regions rather than a universal default, and GFNI whole-codec
    evidence is single-host; both remain explicit requests until a same-binary
    selector screen promotes them.  Within such a fixed,
    offline-calibrated region on a recognized CPU model, AUTO may use another
    startup-qualified immutable table for one operation; unknown models retain
    the baseline.  This changes neither code identity nor wire bytes.
    Explicit backend requests never widen.  A table beyond the process
    baseline is allocated and known-answer-tested on its first explicit
    context request or first eligible AUTO codec setup.  Results are cached.
    For an explicit request, an unavailable ISA returns LEO2_UNSUPPORTED,
    allocation failure returns LEO2_OUT_OF_MEMORY, and a known-answer-test
    failure returns LEO2_INTERNAL_ERROR.  Optional AUTO qualification failure
    retains the reported baseline instead of failing codec creation.
    Independent contexts may select different tables concurrently without
    changing the wire profile or process-global legacy leo_* selection.  NEON
    remains an exact request for the existing active ARM path.  Leopard2 caches
    the first process-wide initialization outcome so concurrent and later
    context constructors report one deterministic result.  An application that
    needs explicit retry after a transient startup failure should call
    leo_init() to completion before its first Leopard2 context construction.
*/
typedef struct leo2_context_options {
    size_t struct_size;
    /* Fixed-width ABI storage; set to one of the leo2_backend values. */
    uint32_t backend;
    uint32_t thread_count; /* zero means library default */
    uint32_t reserved;
} leo2_context_options;

typedef struct leo2_codec_options {
    size_t struct_size;
    uint32_t flags;
    uint32_t reserved;
    /*
        This field was appended in API version 2.  A version-1 struct_size,
        ending immediately before shard_layout, selects NATIVE_V1.
    */
    /* Fixed-width ABI storage; set to one of the leo2_shard_layout values. */
    uint32_t shard_layout;
} leo2_codec_options;

/* Test/diagnostic flags selecting a decoder independently of AUTO dispatch. */
#define LEO2_CODEC_FORCE_GENERIC_DECODE 0x00000001u
#define LEO2_CODEC_FORCE_SPECIALIZED_DECODE 0x00000002u
/* Select the side-sized or retained N-slot specialized kernel for diagnostics. */
#define LEO2_CODEC_FORCE_TILED_DECODE 0x00000004u
#define LEO2_CODEC_FORCE_MATERIALIZED_DECODE 0x00000008u

/*
    These are the fixed V1 batch-item layouts.  A future incompatible extension
    will use a new suffixed type and entry point rather than reading past them.
    A batch may contain at most UINT32_MAX items.
*/
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

/* Accepts any integer so diagnostics remain defined for unknown future codes. */
LEO2_EXPORT const char* leo2_result_string(int result);

/*
    Successful setup returns immutable objects.  A context must outlive all of
    its codecs, and a codec must outlive all of its decode plans.  A non-null
    output handle is set to null before ordinary input validation and remains
    null on failure.  Setup options/presence arrays are immutable and their
    declared spans must not overlap the output handle.  Such overlap returns
    LEO2_OVERLAP without modifying the handle, because clearing it would modify
    immutable input.  An unrepresentable input or output span similarly returns
    LEO2_INVALID_ARGUMENT before the handle is written.  Scalar size-query
    outputs follow the same rule: their complete size_t or uint64_t object span
    must be representable before the function clears or fills it.

    On Linux, a context that explicitly requests AVX2 snapshots the creating
    thread's affinity and the smallest level-three cache visible through that
    mask for cache-sensitive GF16 execution geometry.  AUTO, GFNI, and other
    backends retain the conservative built-in calibration until independently
    qualified.  Missing or inconsistent topology also retains that calibration.
    The snapshot changes only scratch sizing and byte tiling, never the field,
    profile, or wire bytes.  Changing affinity later remains correct, provided
    the later CPUs support the selected ISA, but may be performance-suboptimal.
    Callers that manage affinity should establish it before context creation and
    retain the same mask for ordinary execution and every batch call.  Batch
    workers are created lazily and may grow on a later call, inheriting that
    caller's affinity.
*/
LEO2_EXPORT leo2_result leo2_context_create(
    const leo2_context_options* options,
    leo2_context** context_out);
LEO2_EXPORT void leo2_context_destroy(leo2_context* context);
/*
    Null context introspection returns AUTO or zero.  For an AUTO context,
    leo2_context_backend returns its immutable baseline; deterministic
    per-codec kernel widening described above is intentionally not represented
    as a different context or wire backend.
*/
LEO2_EXPORT leo2_backend leo2_context_backend(const leo2_context* context);
LEO2_EXPORT uint32_t leo2_context_thread_count(const leo2_context* context);

/*
    profile and field intentionally use integer ABI inputs.  Pass one of the
    corresponding leo2_* enum constants.  Keeping untrusted selector values as
    integers until after validation makes negative and unknown future values a
    defined LEO2_INVALID_ARGUMENT result, including under C short-enum modes
    and C++ enum sanitizers.  Introspection returns remain enum-typed because
    the library produces only validated values.
*/
LEO2_EXPORT leo2_result leo2_codec_create(
    leo2_context* context,
    uint32_t original_count,
    uint32_t recovery_count,
    int profile,
    int field,
    const leo2_codec_options* options,
    leo2_codec** codec_out);
LEO2_EXPORT void leo2_codec_destroy(leo2_codec* codec);

/* Null codec introspection returns zero, AUTO, or NATIVE_V1 respectively. */
LEO2_EXPORT uint32_t leo2_codec_original_count(const leo2_codec* codec);
LEO2_EXPORT uint32_t leo2_codec_recovery_count(const leo2_codec* codec);
LEO2_EXPORT uint32_t leo2_codec_parent_count(const leo2_codec* codec);
LEO2_EXPORT uint32_t leo2_codec_padded_side(const leo2_codec* codec);
LEO2_EXPORT leo2_profile leo2_codec_profile(const leo2_codec* codec);
LEO2_EXPORT leo2_field leo2_codec_field(const leo2_codec* codec);
LEO2_EXPORT leo2_shard_layout leo2_codec_shard_layout(const leo2_codec* codec);

/*
    Translate an application payload size to the exact physical shard size.
    NATIVE_V1 returns the payload size subject to the native field granularity.
    GF16_PADDED_ODD_V1 accepts only odd payload sizes and returns B+1.

    Pack/unpack are allocation-free and use memmove semantics, so source and
    destination may overlap.  wire_shard_bytes must equal the queried size.
    An unrepresentable payload or wire span returns LEO2_INVALID_ARGUMENT before
    either range is accessed.  Padded-odd unpack rejects a nonzero systematic
    pad byte.  These helpers are for systematic/original coordinates; parity
    shards must retain all physical bytes and must not be unpacked or truncated.
*/
LEO2_EXPORT leo2_result leo2_codec_wire_shard_bytes(
    const leo2_codec* codec,
    uint64_t payload_bytes,
    uint64_t* wire_shard_bytes_out);
LEO2_EXPORT leo2_result leo2_pack_systematic_shard(
    const leo2_codec* codec,
    uint64_t payload_bytes,
    const void* payload,
    void* wire_shard,
    uint64_t wire_shard_bytes);
LEO2_EXPORT leo2_result leo2_unpack_systematic_shard(
    const leo2_codec* codec,
    uint64_t payload_bytes,
    const void* wire_shard,
    uint64_t wire_shard_bytes,
    void* payload);

/*
    After its complete object span is validated, a scratch-size output is set
    to zero before ordinary argument validation and remains zero on failure.
    Scratch must start at leo2_scratch_alignment() alignment.  Scratch
    and all non-null shard ranges must be mutually disjoint, except that input
    shards may alias other input shards.  Pointer arrays and a batch's item array
    are immutable call metadata: they may share immutable input storage or each
    other, but must not overlap scratch or a writable shard range.  Across a
    batch, immutable input shards may be shared between items; every supplied
    scratch span and writable shard span must be disjoint from every other
    item's readable and writable spans.  The whole batch is preflighted before
    any item executes, so unsupported overlap and ordinary per-item validation
    failures do not partially encode or restore another item.  This
    allocation-free preflight performs O(B^2 * (K+R)) range work for B items;
    the documented UINT32_MAX item-count limit is not silently reduced.
    Unsupported
    metadata overlap is rejected before scratch is modified.  Recovery output
    entries may be null to request a parity subset.  Encode/decode execution
    performs no allocation.  shard_bytes always describes the physical buffers.
    GF8 accepts every positive native length.  GF16 core execution accepts
    positive even physical lengths and returns LEO2_UNSUPPORTED for odd physical
    lengths.  A padded-odd codec also requires the last byte of every systematic
    physical input to be zero.  Use leo2_codec_wire_shard_bytes and the
    pack/unpack helpers to map an odd payload B to physical W=B+1.  A no-loss
    plan remains a zero-scratch no-op: decode and decode-batch do not inspect
    per-item pointer arrays, scratch, shard sizes, or aliases for such a plan.
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
    The ordinary batch entry point above retains its allocation-free
    compatibility preflight, whose work is quadratic in item_count.  These
    additive entry points accept one caller-owned, batch-wide preflight scratch
    span and flatten the same alias contract into an O(M log M) interval sweep,
    where M is proportional to item_count * (K+R).  The queried span is aligned
    to leo2_scratch_alignment(), is independent of shard_bytes and parity-subset
    choices, and may be reused for later calls with the same or smaller item
    count.  Its contents are unspecified on return.

    For fewer than nine items the query returns zero and the execution entry
    point preserves the ordinary compatibility path without inspecting
    preflight_scratch; this deterministic cutoff avoids penalizing small
    batches.  A caller can invoke the ordinary entry point when the query
    returns zero to retain its exact call-wrapper overhead.  Otherwise the full
    supplied preflight scratch span must be
    disjoint from the item array, all pointer arrays, all shard ranges, and
    every per-item scratch span.  An overlap between this span and metadata is
    rejected before this span is modified.  The interval preflight and item
    codec execution allocate no
    memory; a context pool can still start or grow workers lazily as documented
    in leo2_context_options.
*/
LEO2_EXPORT leo2_result leo2_encode_batch_preflight_scratch_size(
    const leo2_codec* codec,
    size_t item_count,
    size_t* scratch_bytes_out);
LEO2_EXPORT leo2_result leo2_encode_batch_with_preflight_scratch(
    const leo2_codec* codec,
    const leo2_encode_batch_item* items,
    size_t item_count,
    void* preflight_scratch,
    size_t preflight_scratch_bytes);

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
/* Null plan introspection reports zero missing originals. */
LEO2_EXPORT uint32_t leo2_decode_plan_missing_original_count(
    const leo2_decode_plan* plan);
LEO2_EXPORT leo2_result leo2_decode_plan_scratch_size(
    const leo2_decode_plan* plan,
    uint64_t shard_bytes,
    size_t* scratch_bytes_out);
/*
    Batch items use this same per-item scratch query.  An internal batch-aware
    kernel choice never requires more per-item storage than the ordinary
    single-item result.  The exact item count affects scheduling and may select
    an equal-or-smaller execution layout; it is not a persistent code or plan
    hint.
*/
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
LEO2_EXPORT leo2_result leo2_decode_plan_batch_preflight_scratch_size(
    const leo2_decode_plan* plan,
    size_t item_count,
    size_t* scratch_bytes_out);
LEO2_EXPORT leo2_result leo2_decode_plan_execute_batch_with_preflight_scratch(
    const leo2_decode_plan* plan,
    const leo2_decode_batch_item* items,
    size_t item_count,
    void* preflight_scratch,
    size_t preflight_scratch_bytes);

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
