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

// In-place two-way LCH butterflies.  x and y must be disjoint shard buffers.
// The multiplier is in Leopard's legacy logarithm representation.  Transform
// callers specialize the kModulus zero-skew sentinel before entering these
// functions, so the backend always receives an ordinary fixed multiplier.
typedef void (*Butterfly2)(
    void* x,
    void* y,
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
    Butterfly2 ff8_ifft_butterfly2;
    Butterfly2 ff8_fft_butterfly2;
    IFFTButterfly2Xor ff8_ifft_butterfly2_xor;
    Butterfly2 ff16_ifft_butterfly2;
    Butterfly2 ff16_fft_butterfly2;
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

const Ops& GetOps();
leo2_backend SelectedBackend();
// Reports the effective public execution backend.  This can differ from the
// fixed-ops table on the existing ARM paths, where legacy native NEON or
// SSE2NEON transform kernels execute around a scalar tail/fallback table.
leo2_backend ExecutionBackend();
bool StartupSelfTestPassed();

}} // namespace leopard::backend
