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

#ifndef CAT_LEOPARD2_DIRECT_TEST_H
#define CAT_LEOPARD2_DIRECT_TEST_H

#include "leopard2.h"

/*
    These hooks are intentionally absent from production builds and the public
    Leopard2 ABI.  Tests set the per-codec mode before concurrent execution.
*/
#ifdef LEO2_ENABLE_TEST_HOOKS

#ifdef __cplusplus
extern "C" {
#endif

typedef enum leo2_test_encode_mode {
    LEO2_TEST_ENCODE_AUTO = 0,
    LEO2_TEST_ENCODE_FORCE_DIRECT = 1,
    LEO2_TEST_ENCODE_FORCE_TRANSFORM = 2
} leo2_test_encode_mode;

/*
    Diagnostic decoder selection used to validate the P=T, N=2P translation
    identity.  A plan captures the codec mode at construction and remains
    immutable if the codec mode is later changed.
*/
typedef enum leo2_test_decode_mode {
    LEO2_TEST_DECODE_AUTO = 0,
    LEO2_TEST_DECODE_FORCE_TRANSLATED_LOW = 1,
    LEO2_TEST_DECODE_FORCE_NATIVE_HIGH = 2
} leo2_test_decode_mode;

/* Exact balanced byte-pass geometry shared by encode and decode scratch.
   Zero aligned bytes is the empty schedule; every nonzero input and requested
   pass size must be a complete 64-byte kernel tile. */
LEO2_EXPORT leo2_result leo2_test_balanced_execution_tiles(
    uint64_t aligned_bytes,
    uint64_t requested_tile_bytes,
    size_t* execution_tile_count_out,
    size_t* maximum_pass_bytes_out);

/* Integer inputs keep deliberately invalid diagnostic probes defined under
   C++ enum sanitizers; valid leo2_test_* enum constants convert implicitly. */
LEO2_EXPORT leo2_result leo2_test_codec_set_encode_mode(
    leo2_codec* codec,
    int mode);

/* Whether immutable binding setup captured the K=1/R=1 copy terminal. */
LEO2_EXPORT int leo2_test_encode_batch_binding_uses_k1_copy(
    const leo2_encode_batch_binding* binding);

LEO2_EXPORT leo2_result leo2_test_codec_set_decode_mode(
    leo2_codec* codec,
    int mode);

LEO2_EXPORT int leo2_test_codec_translated_low_capable(
    const leo2_codec* codec);

LEO2_EXPORT int leo2_test_decode_plan_uses_translated_low(
    const leo2_decode_plan* plan);

/* True when an R=1 plan stores only its missing-original scalar rather than
   owning full original/recovery presence vectors. */
LEO2_EXPORT int leo2_test_decode_plan_is_compact_direct_xor(
    const leo2_decode_plan* plan);

/* Number of immutable source-major coefficient rows prepared for this plan.
   Zero is the production/default-OFF representation. */
LEO2_EXPORT size_t leo2_test_decode_plan_direct_source_rows(
    const leo2_decode_plan* plan);

LEO2_EXPORT int leo2_test_codec_direct_encode_capable(
    const leo2_codec* codec);

LEO2_EXPORT leo2_result leo2_test_codec_encode_path(
    const leo2_codec* codec,
    uint64_t shard_bytes,
    /* Number of non-null requested parity outputs, not a parity prefix. */
    uint32_t requested_recovery_count,
    int* direct_out);

/* Reports the immutable table selected for a transform encode call shape. */
LEO2_EXPORT leo2_result leo2_test_codec_transform_encode_backend(
    const leo2_codec* codec,
    uint64_t shard_bytes,
    uint32_t requested_recovery_count,
    uint32_t requested_recovery_prefix,
    leo2_backend* backend_out);

LEO2_EXPORT void leo2_test_reset_generic_reveal_counts(void);

LEO2_EXPORT uint64_t leo2_test_generic_direct_reveal_shards(void);

LEO2_EXPORT void leo2_test_reset_low_reveal_counts(void);

LEO2_EXPORT uint64_t leo2_test_low_direct_reveal_shards(void);

LEO2_EXPORT uint64_t leo2_test_low_scratch_reveal_shards(void);

LEO2_EXPORT void leo2_test_reset_high_reveal_counts(void);

LEO2_EXPORT uint64_t leo2_test_high_direct_reveal_shards(void);

LEO2_EXPORT uint64_t
leo2_test_high_materialized_direct_reveal_shards(void);

LEO2_EXPORT uint64_t leo2_test_high_scratch_reveal_shards(void);

LEO2_EXPORT void leo2_test_reset_direct_pair_calls(void);

LEO2_EXPORT uint64_t leo2_test_direct_pair_calls(void);

LEO2_EXPORT void leo2_test_reset_direct_four_tiny_calls(void);

LEO2_EXPORT uint64_t leo2_test_direct_four_tiny_calls(void);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* LEO2_ENABLE_TEST_HOOKS */

#ifdef __cplusplus

namespace leopard2_internal {

/*
    Internal structural accounting for immutable decode tables.  This is not
    part of the installed C API: production and hook-library tests use it to
    prove that codec setup retains exactly the metadata reachable through each
    build's dispatch controls.
*/
struct CodecDecodeMetadataInfo
{
    size_t permanent_erased_bytes;
    size_t native_locator_bytes;
    size_t native_factor_bytes;
    size_t translated_permanent_erased_bytes;
    size_t translated_locator_bytes;
    size_t translated_factor_bytes;
};

bool GetCodecDecodeMetadataInfo(
    const leo2_codec* codec,
    CodecDecodeMetadataInfo* info_out);

/*
    Production encode-table and AUTO-route accounting.  This internal query is
    deliberately available without test hooks so an executable linked against
    the ordinary archive can prove whether a candidate selector merely exists
    in a hook build or is actually prepared and selected in production.
*/
struct CodecEncodePathInfo
{
    size_t direct_generator_rows;
    bool auto_direct_selected;
    bool high_t4_batch_binding_selected;
    bool high_t8_vector_selected;
    bool high_t8_partial_binding_selected;
    bool high_t8_two_block_binding_selected;
};

bool GetCodecEncodePathInfo(
    const leo2_codec* codec,
    uint64_t shard_bytes,
    uint32_t requested_recovery_count,
    CodecEncodePathInfo* info_out);

/*
    Deterministic upper byte boundary for the dense legacy-high GF8 AVX2
    T=4 binding shortcut.  Zero means that the public shape is not supported.
    This is an internal policy query, not part of the wire profile or C ABI.
*/
uint64_t HighT4BatchMaximumBytes(
    uint32_t original_count,
    uint32_t recovery_count);

/*
    Setup-only benchmark control for same-executable attribution.  Call this
    immediately after context creation and before constructing any codec.
    This is deliberately internal: It does not alter the public API, codec
    identity, wire profile, or production default.
*/
bool SetContextHighT4BatchBindingEnabledForDiagnostics(
    leo2_context* context,
    bool enabled);

/*
    Process-local benchmark control for the exact-byte K=8/R=3..4 packed
    T=4 terminals.  Change it only while no encode call is executing and
    invoke it before constructing the codec under test.  This remains outside
    the public API and changes no codec or wire identity.
*/
bool SetK8R3R4T4TerminalEnabledForDiagnostics(bool enabled);

/*
    Process-local benchmark control for the exact K=16/R=8/256-byte packed
    terminal.  Invoke before context creation.  This remains outside the
    public API and changes no codec or wire identity.
*/
bool SetK16R8B256TerminalEnabledForDiagnostics(bool enabled);

/*
    Process-local benchmark control for the exact K=9/R=5/256-byte packed
    terminal.  Change it only while no encode call is executing.  This remains
    outside the public API and changes no codec or wire identity.
*/
bool SetK9R5B256TerminalEnabledForDiagnostics(bool enabled);

/*
    Process-local benchmark control for the exact K=9/R=6..8/256-byte packed
    terminals.  Change it only while no encode call is executing.  This
    remains outside the public API and changes no codec or wire identity.
*/
bool SetK9R6R8B256TerminalEnabledForDiagnostics(bool enabled);

/*
    Pattern-dependent structural accounting for an immutable decode plan.
    Empty high-output entries select the mature full transform and are not
    counted as compiled pruned schedules.
*/
struct DecodePlanPrunedScheduleInfo
{
    size_t low_input_plan_count;
    size_t low_output_plan_count;
    size_t high_input_plan_count;
    size_t high_output_plan_count;
};

bool GetDecodePlanPrunedScheduleInfo(
    const leo2_decode_plan* plan,
    DecodePlanPrunedScheduleInfo* info_out);

/* Fail-closed provenance marker for the one-shot no-loss experiment. */
bool OneShotNoLossShortCircuitExperimentEnabled();

/*
    Private benchmark provenance marker.  Candidate/control builds read a
    volatile data initializer rather than compiling different selector code.
*/
bool HighT8OneBlockExtendedEnabled();

/* Provisional marker for the measured one-block extension above 512 bytes. */
bool HighT8OneBlockBeyond512Enabled();

/* Marker for the independently qualified T=8 1024-byte shape extension. */
bool HighT8OneKilobyteExtensionEnabled();

/* Marker for the qualified sub-64-byte T=8 binding. */
bool HighT8TinyBindingEnabled();

/* Marker for the qualified 65..1024-byte ragged T=8 selector. */
bool HighT8RaggedBindingEnabled();

/* Text-layout-neutral marker for equal-rounded GF8/AVX2 multi-loss repair. */
bool EqualRoundedMultiLossEnabled();

/* Equivalent marker for the bounded public one-shot repair path. */
bool OneShotEqualRoundedDirectEnabled();

/* Text-layout-neutral marker for Cauchy setup log reuse. */
bool CauchyLogReuseEnabled();

bool HighT8TwoBlock128192Enabled();

/* Equivalent text-layout-neutral marker for the 320-byte extension. */
bool HighT8TwoBlock320Enabled();

/* Provisional marker for the measured extension above 320 bytes. */
bool HighT8TwoBlockExtendedEnabled();

#ifdef LEO2_ENABLE_TEST_HOOKS
/*
    Test-only storage accounting for presence-dependent no-op plan state.
    Capacities distinguish an empty logical view from retained allocations.
*/
bool GetDecodePlanPresenceStorageInfo(
    const leo2_decode_plan* plan,
    size_t* original_capacity_out,
    size_t* recovery_capacity_out,
    size_t* erased_capacity_out);
#endif

} // namespace leopard2_internal

#endif /* __cplusplus */

#endif /* CAT_LEOPARD2_DIRECT_TEST_H */
