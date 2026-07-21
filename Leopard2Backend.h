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

#pragma once

#include "leopard2.h"

#include <stdint.h>

namespace leopard { namespace backend {

typedef uint8_t (*FF8MultiplyLog)(uint8_t value, uint8_t multiplier_log);
typedef uint16_t (*FF16MultiplyLog)(uint16_t value, uint16_t multiplier_log);

struct InitializeArgs
{
    FF8MultiplyLog ff8_multiply_log;
    FF16MultiplyLog ff16_multiply_log;
};

typedef void (*FixedMultiply)(
    void* destination,
    const void* source,
    uint16_t multiplier_log,
    uint64_t byte_count);

typedef void (*XorMemory)(
    void* destination,
    const void* source,
    uint64_t byte_count);

// One-pass accumulation of two read-only inputs:
//
//     destination[i] ^= source0[i] ^ source1[i]
//
// The destination range must be disjoint from both source ranges.  The two
// read-only sources may alias one another, matching the public API's allowed
// input-input aliasing.  A zero byte count performs no access.
typedef void (*XorMemory2To1)(
    void* destination,
    const void* source0,
    const void* source1,
    uint64_t byte_count);

// Copy initial_source to destination, then XOR every non-null entry in
// sources[0..source_count).  The destination must be disjoint from every
// source range; read-only source ranges may alias one another.  This coarse
// operation keeps the R=1 reduction loop inside the selected ISA translation
// unit instead of paying one indirect Ops dispatch per source pair.
typedef void (*XorMemorySources)(
    void* destination,
    const void* initial_source,
    const void* const* sources,
    uint32_t source_count,
    uint64_t byte_count);

// Four independent in-place XOR pairs.  Destination ranges must be pairwise
// disjoint and disjoint from all source ranges.  Read-only source ranges may
// alias one another, matching the public API's allowed input-input aliasing.
// Keeping this grouped operation in Ops prevents VectorXOR from silently
// falling back to the process-default ISA for an explicitly lower context.
typedef void (*XorMemory4)(
    void* destination0,
    const void* source0,
    void* destination1,
    const void* source1,
    void* destination2,
    const void* source2,
    void* destination3,
    const void* source3,
    uint64_t byte_count);

// In-place two-way LCH butterflies.  x and y must be disjoint shard buffers.
// The multiplier is in Leopard's legacy logarithm representation.  Transform
// callers specialize the kModulus zero-skew sentinel before entering these
// functions, so the backend always receives an ordinary fixed multiplier.
typedef void (*Butterfly2)(
    void* x,
    void* y,
    uint16_t multiplier_log,
    uint64_t byte_count);

// Out-of-place forward two-way LCH butterfly.  The input ranges are read-only
// and must be disjoint from the two disjoint output ranges.  The field's
// all-ones logarithm value is accepted as the zero-skew sentinel: in that case
// the multiplication is suppressed but the XOR edge is retained.  This form
// lets low-rate encoding preserve its coefficient block while combining the
// former coefficient copy with the first transform layer.
typedef void (*FFTButterfly2Out)(
    const void* x_input,
    const void* y_input,
    void* x_output,
    void* y_output,
    uint16_t multiplier_log,
    uint64_t byte_count);

// Accumulating inverse two-way butterfly used by the encoders.  All four
// buffers must be pairwise disjoint.  Inputs are read-only; outputs are XORed
// with the inverse-butterfly result.  Transform callers specialize the
// field's all-ones zero-skew sentinel before entering this function.
typedef void (*IFFTButterfly2Xor)(
    const void* x_input,
    const void* y_input,
    void* x_output,
    void* y_output,
    uint16_t multiplier_log,
    uint64_t byte_count);

// In-place fused two-layer LCH butterfly over four disjoint shard buffers.
// The three multipliers correspond to pairs (0,1), (2,3), and (0,2)/(1,3).
// UINT8_MAX and UINT16_MAX are the GF8/GF16 zero-skew sentinels respectively;
// a sentinel suppresses the fixed multiplication but retains the XOR edge.
typedef void (*Butterfly4)(
    void* value0,
    void* value1,
    void* value2,
    void* value3,
    uint16_t multiplier_log01,
    uint16_t multiplier_log23,
    uint16_t multiplier_log02,
    uint64_t byte_count);

// Executes one complete contiguous radix-four transform group.  The field
// scheduler has already proved that work[0..4*distance) names disjoint shard
// buffers; each lane i applies the same three multipliers to
// { i, i + distance, i + 2*distance, i + 3*distance }.  Keeping this loop in
// the ISA-specific translation unit reduces one indirect Ops dispatch per
// butterfly to one dispatch per group while preserving per-context backend
// selection.  prefer_fused records a caller-selected, measured whole-stage
// policy.  When false, each backend retains its mature per-butterfly size
// policy.
typedef void (*Butterfly4Range)(
    void* const* work,
    unsigned distance,
    uint16_t multiplier_log01,
    uint16_t multiplier_log23,
    uint16_t multiplier_log02,
    uint64_t byte_count,
    bool prefer_fused);

// GF8 encoder variant that accumulates a completed inverse radix-four group
// into four corresponding output ranges.  The output ranges are XORed and
// must be disjoint from the work ranges.  Work contents are unspecified after
// the call: a backend may transform them in place or consume them read-only
// while fusing the final butterfly with output accumulation.
typedef void (*IFFTButterfly4XorRange)(
    void* const* work,
    void* const* xor_output,
    unsigned distance,
    uint16_t multiplier_log01,
    uint16_t multiplier_log23,
    uint16_t multiplier_log02,
    uint64_t byte_count);

// Out-of-place forward fused two-layer LCH butterfly.  All input ranges are
// read-only, all output ranges are disjoint, and inputs must not overlap any
// output.  Sentinel behavior matches Butterfly4.  Sources are loaded directly
// into registers and transformed before the first store; this is not a
// copy-then-in-place compatibility wrapper.
typedef void (*FFTButterfly4Out)(
    const void* input0,
    const void* input1,
    const void* input2,
    const void* input3,
    void* output0,
    void* output1,
    void* output2,
    void* output3,
    uint16_t multiplier_log01,
    uint16_t multiplier_log23,
    uint16_t multiplier_log02,
    uint64_t byte_count);

// Out-of-place inverse companion.  The storage and sentinel contract is the
// same as FFTButterfly4Out; only the LCH butterfly direction differs.  High-
// rate encoding uses this to consume four caller shards directly into the
// first two inverse layers instead of copying them through a full workspace
// pass first.
typedef FFTButterfly4Out IFFTButterfly4Out;

// Complete one-block legacy-high GF8 transform.  The caller has already
// established that data[0..data_count) and work[0..side) are disjoint and that
// every transmitted parity coordinate is requested.  side_and_flags stores
// side in the low bits and sets kFF8HighEncodeShortenedInput when data_count is
// side - 1; otherwise data_count is side.  Keeping this one-bit fact in the
// existing size argument preserves the six-register hot-call shape on x86-64.
// A missing final data coordinate is the shortened known zero;
// a punctured final parity coordinate is still evaluated into the caller's
// scratch-backed work row.  inverse_skew names the message-coset skew table,
// while forward_skew names the parity-coset table.  Keeping the regular
// transform traversal in an ISA translation unit removes one indirect Ops
// dispatch per butterfly range without changing the wire profile or arithmetic
// order.  This callback is optional.
static const uint32_t kFF8HighEncodeShortenedInput = 0x80000000U;

typedef void (*FF8HighEncodeOneBlock)(
    const void* const* data,
    void* const* work,
    uint32_t side_and_flags,
    const uint8_t* inverse_skew,
    const uint8_t* forward_skew,
    uint64_t byte_count);

// Complete dense legacy-high GF8 encode for the two- and four-coordinate
// redundancy transforms.  The implementation consumes all original shards
// directly, accumulates each shifted inverse transform into work[0..side),
// and performs the final parity-coset transform in place.  This avoids the
// compatibility path's per-message-block copy through work[side..2*side).
//
// side must be 2 or 4.  data[0..original_count) are read-only and disjoint
// from the pairwise-disjoint work[0..side) outputs.  inverse_skew points at
// the first message coset (the field table's +side entry), while forward_skew
// points at the parity coset.  The callback is optional and may be used for a
// contiguous requested parity prefix: it materializes all side parent-code
// parity coordinates, placing unrequested or punctured suffix coordinates in
// caller-provided scratch slots.
typedef void (*FF8HighEncodeSmall)(
    const void* const* data,
    uint32_t original_count,
    void* const* work,
    uint32_t side,
    const uint8_t* inverse_skew,
    const uint8_t* forward_skew,
    uint64_t byte_count);

// GF8 AVX2 boundary operation.  Applies four independent nonzero fixed
// multipliers to the selected input rows and immediately executes the first
// two inverse LCH layers.  weight_log values are ordinary GF8 logarithms: both
// zero and 255 denote the multiplicative identity.  This differs deliberately
// from the butterfly logs, where 255 remains the zero-skew sentinel.  A
// cleared bit in live_mask supplies mathematical zero and permits the
// corresponding input pointer to be null.
//
// Outputs must be pairwise disjoint.  The inputs may either be wholly disjoint
// from every output (and may alias one another), or every output must exactly
// equal its corresponding input for a four-row in-place operation.  Partial
// and cross-row overlap are not supported.  Implementations load all live
// inputs for one byte/symbol group before storing any result.
typedef void (*WeightedIFFTButterfly4)(
    const void* input0,
    const void* input1,
    const void* input2,
    const void* input3,
    void* output0,
    void* output1,
    void* output2,
    void* output3,
    uint16_t weight_log0,
    uint16_t weight_log1,
    uint16_t weight_log2,
    uint16_t weight_log3,
    uint8_t live_mask,
    uint16_t multiplier_log01,
    uint16_t multiplier_log23,
    uint16_t multiplier_log02,
    uint64_t byte_count);

// Shared debug/test oracle for the weighted boundary's deliberately narrow
// alias contract.  A zero-byte operation accesses no range.
bool IsWeightedIFFTButterfly4AliasingValid(
    const void* const inputs[4],
    const void* const outputs[4],
    uint8_t live_mask,
    uint64_t byte_count);

// This table is private to the implementation and immutable.  A backend owns
// any tables referenced by its functions and publishes this object only after
// initialization and the startup known-answer tests have succeeded.
struct Ops
{
    leo2_backend kind;
    const char* name;
    FixedMultiply ff8_multiply;
    FixedMultiply ff8_multiply_add;
    FixedMultiply ff16_multiply;
    FixedMultiply ff16_multiply_add;
    XorMemory xor_memory;
    XorMemory2To1 xor_memory_2to1;
    XorMemorySources xor_memory_sources;
    XorMemory4 xor_memory4;
    Butterfly2 ff8_ifft_butterfly2;
    Butterfly2 ff8_fft_butterfly2;
    FFTButterfly2Out ff8_fft_butterfly2_out;
    IFFTButterfly2Xor ff8_ifft_butterfly2_xor;
    Butterfly4 ff8_ifft_butterfly4;
    Butterfly4 ff8_fft_butterfly4;
    IFFTButterfly4Out ff8_ifft_butterfly4_out;
    // Optional: the qualified AVX2-family backends currently supply this.
    WeightedIFFTButterfly4 ff8_weighted_ifft_butterfly4;
    FFTButterfly4Out ff8_fft_butterfly4_out;
    Butterfly4Range ff8_ifft_butterfly4_range;
    Butterfly4Range ff8_fft_butterfly4_range;
    IFFTButterfly4XorRange ff8_ifft_butterfly4_xor_range;
    Butterfly2 ff16_ifft_butterfly2;
    Butterfly2 ff16_fft_butterfly2;
    FFTButterfly2Out ff16_fft_butterfly2_out;
    IFFTButterfly2Xor ff16_ifft_butterfly2_xor;
    Butterfly4 ff16_ifft_butterfly4;
    Butterfly4 ff16_fft_butterfly4;
    IFFTButterfly4Out ff16_ifft_butterfly4_out;
    FFTButterfly4Out ff16_fft_butterfly4_out;
    Butterfly4Range ff16_ifft_butterfly4_range;
    Butterfly4Range ff16_fft_butterfly4_range;
    FF8HighEncodeOneBlock ff8_high_encode_one_block;
    FF8HighEncodeSmall ff8_high_encode_small;
};

struct X86Features
{
    bool ssse3;
    bool avx2;
    bool avx512;
};

struct X86ProcessorIdentity
{
    bool authentic_amd;
    uint16_t family;
    uint16_t model;
};

X86ProcessorIdentity ClassifyX86Processor(
    uint32_t vendor_ebx,
    uint32_t vendor_edx,
    uint32_t vendor_ecx,
    uint32_t leaf1_eax);
X86ProcessorIdentity DetectX86Processor();
bool IsCalibratedAutoAVX512EncodeProcessor(
    const X86ProcessorIdentity& identity);
bool IsCalibratedAutoAVX512EncodeHost();

// Pure feature classifier used by both the runtime probe and synthetic tests.
X86Features ClassifyX86Features(
    uint32_t maximum_basic_leaf,
    uint32_t leaf1_ecx,
    uint32_t leaf7_ebx,
    uint64_t xcr0);
X86Features DetectX86Features();

const Ops* InitializeScalar(const InitializeArgs& args);

#if defined(LEO2_HAVE_SSSE3_BACKEND)
const Ops* InitializeSSSE3(const InitializeArgs& args);
#endif

#if defined(LEO2_HAVE_AVX2_BACKEND)
const Ops* InitializeAVX2(const InitializeArgs& args);
// Immutable nibble-table storage shared by the AVX2 and AVX-512VL codegen
// variants.  The erased pointer types keep the private table layouts local to
// their implementation translation unit.
const void* GetAVX2FF8Tables();
const void* GetAVX2FF16Tables();
#endif

#if defined(LEO2_HAVE_AVX512_BACKEND)
const Ops* InitializeAVX512(const InitializeArgs& args);
#endif

// Called once after both legacy fields have initialized their logarithm tables.
// Returns false on allocation or known-answer-test failure.
bool Initialize(const InitializeArgs& args);

// Returns the immutable process default.  Byte-heavy Leopard2 entry points
// pass their context table explicitly instead of consulting thread-local state.
const Ops& GetOps();
const Ops& GetDefaultOps();
enum QualificationStatus
{
    QualificationAvailable = 0,
    QualificationUnavailable,
    QualificationOutOfMemory,
    QualificationSelfTestFailed
};

// AUTO resolves to the process default.  On first use, an explicit lower value
// initializes and known-answer-tests its immutable table under a setup lock.
// A failed or unavailable request is cached and returns null deterministically.
const Ops* GetQualifiedOps(
    leo2_backend requested,
    QualificationStatus* status = NULL);
leo2_backend SelectedBackend();
// Reports the effective public execution backend.  This can differ from the
// fixed-ops table on the existing ARM paths, where legacy native NEON or
// SSE2NEON transform kernels execute around a scalar tail/fallback table.
leo2_backend ExecutionBackend();
bool StartupSelfTestPassed();
// Distinguishes backend allocation and known-answer-test failures for the new
// API while the legacy leo_init() contract continues to report Platform.
QualificationStatus StartupQualificationFailure();

#ifdef LEO2_ENABLE_TEST_HOOKS
// Deterministic one-shot setup faults.  These hooks are compiled only into
// test-enabled archives; production archives contain neither the selector nor
// any branches in setup or byte-execution paths.
enum TestSetupFault
{
    TestSetupFaultNone = 0,
    TestSetupFaultScalarFF8Allocation,
    TestSetupFaultScalarFF16Allocation,
    TestSetupFaultSSSE3FF8Allocation,
    TestSetupFaultSSSE3FF16Allocation,
    TestSetupFaultAVX2FF8Allocation,
    TestSetupFaultAVX2FF16Allocation,
    TestSetupFaultScalarKAT,
    TestSetupFaultSSSE3KAT,
    TestSetupFaultAVX2KAT,
    TestSetupFaultAVX512FF8Allocation,
    TestSetupFaultAVX512FF16Allocation,
    TestSetupFaultAVX512KAT
};

struct TestBackendState
{
    bool ff8_published;
    bool ff16_published;
    bool qualified;
    QualificationStatus failure;
    uint64_t ff8_bytes;
    uint64_t ff16_bytes;
};

void TestSetSetupFault(TestSetupFault fault);
bool TestSetupFaultPending();
unsigned TestSetupFaultConsumptions();
bool TestShouldFailAllocation(leo2_backend backend, bool ff16);
leo2_backend TestDefaultBackendForHost();
bool TestBackendCanQualifyForHost(leo2_backend backend);
bool TestGetBackendState(leo2_backend backend, TestBackendState* state);

void TestGetScalarTableState(TestBackendState* state);
# if defined(LEO2_HAVE_SSSE3_BACKEND)
void TestGetSSSE3TableState(TestBackendState* state);
# endif
# if defined(LEO2_HAVE_AVX2_BACKEND)
void TestGetAVX2TableState(TestBackendState* state);
# endif
# if defined(LEO2_HAVE_AVX512_BACKEND)
void TestGetAVX512TableState(TestBackendState* state);
# endif

// Test-only injection point for a copied/tracing immutable ops table.  The
// caller must restore the qualified table before destroying the tracing table.
void TestSetContextOps(leo2_context* context, const Ops* ops);
#endif

}} // namespace leopard::backend
