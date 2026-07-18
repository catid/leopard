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

// Accumulating inverse two-way butterfly used by the GF8 encoder.  All four
// buffers must be pairwise disjoint.  Inputs are read-only; outputs are XORed
// with the inverse-butterfly result.
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
    XorMemory4 xor_memory4;
    Butterfly2 ff8_ifft_butterfly2;
    Butterfly2 ff8_fft_butterfly2;
    FFTButterfly2Out ff8_fft_butterfly2_out;
    IFFTButterfly2Xor ff8_ifft_butterfly2_xor;
    Butterfly4 ff8_ifft_butterfly4;
    Butterfly4 ff8_fft_butterfly4;
    FFTButterfly4Out ff8_fft_butterfly4_out;
    Butterfly2 ff16_ifft_butterfly2;
    Butterfly2 ff16_fft_butterfly2;
    FFTButterfly2Out ff16_fft_butterfly2_out;
    Butterfly4 ff16_ifft_butterfly4;
    Butterfly4 ff16_fft_butterfly4;
    FFTButterfly4Out ff16_fft_butterfly4_out;
};

struct X86Features
{
    bool ssse3;
    bool avx2;
};

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
    TestSetupFaultAVX2KAT
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
bool TestGetBackendState(leo2_backend backend, TestBackendState* state);

void TestGetScalarTableState(TestBackendState* state);
# if defined(LEO2_HAVE_SSSE3_BACKEND)
void TestGetSSSE3TableState(TestBackendState* state);
# endif
# if defined(LEO2_HAVE_AVX2_BACKEND)
void TestGetAVX2TableState(TestBackendState* state);
# endif

// Test-only injection point for a copied/tracing immutable ops table.  The
// caller must restore the qualified table before destroying the tracing table.
void TestSetContextOps(leo2_context* context, const Ops* ops);
#endif

}} // namespace leopard::backend
