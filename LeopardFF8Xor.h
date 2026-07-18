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

#include "LeopardCommon.h"

#ifdef LEO_HAS_FF8

namespace leopard { namespace ff8xor {


//------------------------------------------------------------------------------
// Field constants

typedef uint8_t ffe_t;

static const unsigned kBits = 8;
static const unsigned kOrder = 256;
static const ffe_t kModulus = 255;
static const unsigned kPolynomial = 0x11D;


//------------------------------------------------------------------------------
// Kernel inspection and test controls

enum class KernelMode
{
    Auto,
    Portable,
    Simd128,
    Avx2,
    // AVX-512VL uses 256-bit YMM values while exposing the extended register
    // file.  Avx512Zmm uses 512-bit values.  Both are forced experiment modes;
    // Auto remains on the established AVX2/128/portable selection until a
    // complete-codec benchmark justifies changing it.
    Avx512VL,
    Avx512Zmm
};

// The four-buffer path is intentionally independent of KernelMode.  Auto
// continues to resolve to the established backend; these controls select the
// lowering only when a caller explicitly selects AVX-512 ZMM and a complete
// two-layer radix-4 unit is safe to fuse.
enum class FourBufferMode
{
    Disabled,
    Xor2,
    Xor3
};

struct CircuitStatistics
{
    unsigned MinimumGateCount;
    unsigned MaximumGateCount;
    double AverageGateCount;
    unsigned MinimumDepth;
    unsigned MaximumDepth;
    double AverageDepth;
};

struct CircuitCost
{
    unsigned GateCount;
    unsigned Depth;
};

// Statistics for the most recent transform on this thread.  Payload byte
// counts are deterministic scheduled traffic estimates (loads plus stores),
// not hardware performance-counter measurements.
struct MaterializationStatistics
{
    uint64_t DeferredZeroFills;
    uint64_t AddedZeroFills;
    uint64_t ButterfliesSkipped;
    uint64_t ButterfliesReduced;
    uint64_t XorsSkipped;
    uint64_t XorsReplacedByCopies;
    uint64_t IdentityOperationsElided;
    int64_t EstimatedPayloadBytesElided;
};

struct FourBufferStatistics
{
    uint64_t FusedUnits;
    uint64_t EstimatedPayloadBytesElided;
};

void SetKernelMode(KernelMode mode);
KernelMode GetKernelMode();
KernelMode GetActiveKernelMode();
bool IsKernelModeAvailable(KernelMode mode);
const char* GetKernelBackendName();
void SetFourBufferMode(FourBufferMode mode);
FourBufferMode GetFourBufferMode();

// Pure register-level predicate used by the runtime detector and tests.  Both
// the processor feature bits and the OS-enabled XMM/YMM/ZMM/opmask state are
// required before entering either AVX-512-targeted width.
bool IsAVX512Supported(
    uint32_t maximum_basic_leaf,
    uint32_t leaf1_ecx,
    uint32_t leaf7_ebx,
    uint64_t xcr0);

const char* GetCircuitChecksum();
const char* GetCircuitCostProfileId();
const char* GetCircuitCostProfileChecksum();
CircuitStatistics GetMultiplyCircuitStatistics();
CircuitStatistics GetFFTCircuitStatistics();
CircuitStatistics GetIFFTCircuitStatistics();
CircuitCost GetMultiplyCircuitCost(ffe_t log_multiplier);
CircuitCost GetFFTCircuitCost(ffe_t skew);
CircuitCost GetIFFTCircuitCost(ffe_t skew);
MaterializationStatistics GetLastMaterializationStatistics();
void ResetMaterializationStatistics();
FourBufferStatistics GetLastFourBufferStatistics();
void ResetFourBufferStatistics();

// Whole-buffer field addition.  This follows the selected experimental kernel
// mode and is also used by the native API's M=1 and formal-derivative paths so
// forced Portable never falls through to a -march=native packed helper.
void XorBuffer(
    uint64_t buffer_bytes,
    void* destination,
    const void* source);

// Internal/test-visible batched addition used by encoder accumulation and the
// formal derivative.  SIMD mode resolution occurs once for the full batch.
// A source may equal its destination within one pair.  Different pairs are
// independent and their memory ranges must not overlap.
void XorBuffers(
    uint64_t buffer_bytes,
    unsigned count,
    void** destination,
    void** source);

// The buffer byte count includes all eight equal contiguous planes.
// Multiplier logarithms 0 and 255 both mean multiplication by one.  Exact
// in-place identity performs no payload access; out-of-place identity has
// memmove-compatible overlap semantics.
void MultiplyBuffer(
    uint64_t buffer_bytes,
    void* destination,
    const void* source,
    ffe_t log_multiplier);

// Unlike multiplier logarithm 255, butterfly skew 255 is a sentinel that
// omits the multiply-add term: x is unchanged and y is XORed with x.
void FFTButterflyBuffer(
    uint64_t buffer_bytes,
    void* x,
    void* y,
    ffe_t skew);

void IFFTButterflyBuffer(
    uint64_t buffer_bytes,
    void* x,
    void* y,
    ffe_t skew);

unsigned GetLastLocatorShiftForTesting();
// -1 restores automatic gate/depth minimization; 0..254 forces a common shift.
void SetLocatorShiftForTesting(int shift);

// Pure exact selector hook.  Logarithm 255 is canonicalized to multiply-log
// zero; inverse entries score -(log + shift).  The objective is minimum total
// gates, then minimum total depth, then lowest numeric shift.  A zero count
// permits null pointers.  Invalid pointers/count return the sentinel 255.
unsigned SelectLocatorShiftForTesting(
    const ffe_t* logarithms,
    const bool* inverse,
    unsigned count);

// Retained shift-major implementation of the same objective for reproducible
// selector microbenchmarks.  Tests compare both hooks with a separate oracle.
unsigned SelectLocatorShiftReferenceForTesting(
    const ffe_t* logarithms,
    const bool* inverse,
    unsigned count);

// Exercise the logical-state butterfly planner directly.  The caller must
// provide payloads matching the declared zero/equality relation.  Deferred
// zero outputs are materialized before return for byte-exact comparison.
void TrackedButterflyBufferForTesting(
    uint64_t buffer_bytes,
    void* x,
    void* y,
    bool inverse,
    bool x_is_zero,
    bool y_is_zero,
    bool equal_nonzero,
    ffe_t skew);

// Compare the direct combined boundary against an independently staged formal
// derivative plus top sentinel FFT.  count must be a power of two in [2, 256]
// and the buffers must be distinct; the operation is naturally in place.
void FormalDerivativeTopFFTForTesting(
    uint64_t buffer_bytes,
    unsigned count,
    void** work);

// Apply one generated four-buffer map through the same runtime gates used by
// production transforms.  False reports an unavailable mode, unknown tuple,
// or a plane size requiring the narrower two-way fallback.
bool FourBufferButterflyBufferForTesting(
    uint64_t buffer_bytes,
    void* a,
    void* b,
    void* c,
    void* d,
    bool inverse,
    ffe_t skew01,
    ffe_t skew23,
    ffe_t skew02);


//------------------------------------------------------------------------------
// Backend API

bool Initialize();
bool IsInitialized();

void ReedSolomonEncode(
    uint64_t buffer_bytes,
    unsigned original_count,
    unsigned recovery_count,
    unsigned m,
    const void* const* data,
    void** work);

void ReedSolomonDecode(
    uint64_t buffer_bytes,
    unsigned original_count,
    unsigned recovery_count,
    unsigned m,
    unsigned n,
    const void* const* original,
    const void* const* recovery,
    void** work);


}} // namespace leopard::ff8xor

#endif // LEO_HAS_FF8
