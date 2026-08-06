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

#ifndef LEO2_EXPERIMENT_HIGH_T16_Q2_B64_FUSED
#define LEO2_EXPERIMENT_HIGH_T16_Q2_B64_FUSED 0
#endif

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

#if LEO2_EXPERIMENT_LOW_P32_B64_TERMINAL
LEO2_EXPORT void leo2_test_reset_low_p32_b64_terminal_calls(void);

LEO2_EXPORT uint64_t leo2_test_low_p32_b64_terminal_calls(void);

LEO2_EXPORT void leo2_test_reset_low_p128_b64_terminal_calls(void);

LEO2_EXPORT uint64_t leo2_test_low_p128_b64_terminal_calls(void);
#endif

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
    size_t translated_full_loss_locator_bytes;
    size_t translated_factor_bytes;
    size_t codec_direct_repair_generator_bytes;
    size_t codec_k8r4_terminal_cache_bytes;
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
    Context-local setup-only attribution control for the pure-AVX2 dense GF8
    Walsh locator.  Call immediately after context creation and before codec
    or plan construction.  It changes neither public ABI nor wire identity.
*/
bool SetContextGF8AVX2WalshLocatorEnabledForDiagnostics(
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

/* Process-local same-executable attribution for the packed K=9..16,
   R=5..8, B64 two-block T=8 terminal. */
bool SetHighT8TwoBlockB64PackedTerminalEnabledForDiagnostics(bool enabled);

/* Process-local same-executable attribution for the packed B256 T=8
   tail/two-block terminal: K=10/R=5..8, K=11/R=5, and K=13..16/R=5..8,
   excluding K=16/R=8. */
bool SetHighT8TwoBlockB256PackedTerminalEnabledForDiagnostics(bool enabled);

/*
    Process-local same-executable control for the dense K=R=5..8 T=8
    packed-slab terminal.  Change it only while no encode call is executing.
    This is benchmark instrumentation, not public API or wire identity.
*/
bool SetT8FullParityTerminalEnabledForDiagnostics(bool enabled);

#if LEO2_EXPERIMENT_HIGH_T16_Q2_B64_FUSED
/*
    Same-executable benchmark control for the default-on
    K=17..32/R=9..16/1..64-byte packed AVX2 terminal. Change it only while no
    encode call is executing.  It changes neither public ABI nor wire identity.
*/
bool SetHighT16Q2B64FusedEnabledForDiagnostics(bool enabled);
#endif

/*
    Process-local benchmark control for the exact balanced 64-byte packed
    terminals.  Change it only while no encode call is executing.  This
    remains outside the public API and changes no codec or wire identity.
*/
bool SetBalancedB64TerminalEnabledForDiagnostics(bool enabled);

bool SetHighT16PreparedTerminalEnabledForDiagnostics(bool enabled);

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
    Same-executable R=1 small-reduction attribution.  Mode zero preserves the
    production policy and mode one enables only the bounded GF8/AVX2 candidate
    cells.  A codec snapshots the process-local mode during construction, so
    changing it cannot alter an existing immutable codec.  Serialize each
    set-and-create pair: concurrent codec construction can otherwise observe
    another benchmark lane's mode.  Values above one are rejected.  This is
    an internal benchmark control, not public ABI.
*/
bool SetR1SmallReductionModeForDiagnostics(unsigned mode);

/*
    Select the fixed-size AVX2 R=1 reducer for codecs created after this call.
    Existing codecs retain their immutable setup-time classification.
    This is diagnostic process-global state; callers must serialize changes
    with codec creation when deterministic attribution is required.
*/
bool SetR1FixedAVX2XorModeForDiagnostics(unsigned mode);

/* Build provenance for the isolated fixed-size AVX2 R=1 candidate. */
bool R1FixedAVX2XorCandidateEnabledForDiagnostics();

/*
    Same-executable attribution for ephemeral one-shot decode-plan setup.
    Mode zero prepares the reusable AUTO metadata superset; mode one uses the
    call's exact shard byte count to prepare only the transform route that
    execution selects; mode two additionally omits
    high-rate pruned schedules in the bounded tiny-GF8/AVX2 attribution
    region; mode three (the production default) extends that policy to the
    translated low-rate execution view.  Public reusable plans are unaffected.
    Change this only while no one-shot decode is executing.
*/
bool SetOneShotPlanSetupModeForDiagnostics(unsigned mode);

/*
    Same-executable attribution for the exact GF8/AVX2 Algorithm 4
    P=32/N=64/B=64 terminal.  Enabled/disabled arm a one-call route probe,
    reported as normalized mode words one/two; finishing the probe restores
    the corresponding timing state and removes route-accounting overhead.
    Both modes execute from identical text.  Change this only while no decode
    is executing; this is an internal benchmark control, not public ABI.
*/
bool SetLowP32B64TerminalEnabledForDiagnostics(bool enabled);

unsigned LowP32B64TerminalModeForDiagnostics();

// Read the calling thread's armed route probe, then disarm it before timing.
bool LowP32B64TerminalRouteSelectedForDiagnostics();

bool FinishLowP32B64TerminalRouteProbeForDiagnostics();

/*
    Same-executable attribution for the bounded GF8/AVX2 Algorithm 4
    P=128/N=256/B=64 terminal.  Its exact target cells are production-enabled
    after qualification; the diagnostic state and route probe remain
    independent of the P32 terminal.
*/
bool SetLowP128B64TerminalEnabledForDiagnostics(bool enabled);

unsigned LowP128B64TerminalModeForDiagnostics();

bool LowP128B64TerminalRouteSelectedForDiagnostics();

bool FinishLowP128B64TerminalRouteProbeForDiagnostics();

/*
    Same-executable attribution for direct final-layer output from dense
    partial GF8/AVX2 P=16 LOW_V1 parity blocks.  Enabled/disabled arm a
    one-call route probe; finishing it normalizes the mode before timing.
*/
bool SetLowP16PartialDirectOutputEnabledForDiagnostics(bool enabled);

unsigned LowP16PartialDirectOutputModeForDiagnostics();

bool LowP16PartialDirectOutputRouteSelectedForDiagnostics();

bool FinishLowP16PartialDirectOutputRouteProbeForDiagnostics();

/*
    Construct the ephemeral transform plan for the current diagnostic mode
    and exact shard byte count.  The factory includes the same no-loss prefix
    probe used by the public wrapper so it does not duplicate that prefix in
    plan validation.  Terminal and direct-repair shapes return
    LEO2_UNSUPPORTED.  This exists only so benchmarks can report transient
    pattern setup separately from execution.  The plan may execute only at
    the byte count supplied here.  Destroy it with leo2_decode_plan_destroy().
*/
leo2_result CreateOneShotTransformPlanForDiagnostics(
    const leo2_codec* codec,
    uint64_t shard_bytes,
    const uint8_t* original_present,
    const uint8_t* recovery_present,
    leo2_decode_plan** plan_out);

leo2_result ExecuteOneShotTransformPlanForDiagnostics(
    const leo2_decode_plan* plan,
    uint64_t shard_bytes,
    const uint8_t* original_present,
    const uint8_t* recovery_present,
    const void* const* original,
    const void* const* recovery,
    void* const* restored_original,
    void* scratch,
    size_t scratch_bytes);

enum R1ReductionPath
{
    kR1ReductionNotApplicable = 0,
    kR1ReductionK1Copy,
    kR1ReductionK2Terminal,
    kR1ReductionPairwise,
    kR1ReductionDense,
    kR1ReductionCoarse,
    kR1ReductionFusedFinal,
    kR1ReductionGroup4,
    kR1ReductionFixedAVX2
};

struct CodecR1ReductionPathInfo
{
    bool small_reduction_mode_enabled;
    R1ReductionPath encode_path;
    R1ReductionPath decode_path;
};

/*
    Reports the reduction body selected by encode execution and by a one-loss
    direct-XOR decode plan for this byte count.  A no-loss or non-direct
    decode plan may not execute the reported decode reduction.
*/
bool GetCodecR1ReductionPathInfo(
    const leo2_codec* codec,
    uint64_t shard_bytes,
    CodecR1ReductionPathInfo* info_out);

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

/*
    Test-only capacity accounting for the incremental direct-repair view held
    by a reusable dual plan.  The byte result excludes allocator bookkeeping
    and the fixed-size opaque plan object, which are shared by both controls.
*/
bool GetDecodePlanDirectStorageInfo(
    const leo2_decode_plan* plan,
    size_t* retained_bytes_out,
    size_t* term_count_out,
    size_t* source_row_count_out);

/* True only when the direct rows were derived from the already-required
   transform locator rather than the independent matrix fallback. */
bool DecodePlanUsesLocatorDirectTerms(const leo2_decode_plan* plan);
#endif

} // namespace leopard2_internal

#endif /* __cplusplus */

#endif /* CAT_LEOPARD2_DIRECT_TEST_H */
