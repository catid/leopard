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
    Avx2
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

void SetKernelMode(KernelMode mode);
KernelMode GetKernelMode();
KernelMode GetActiveKernelMode();
bool IsKernelModeAvailable(KernelMode mode);
const char* GetKernelBackendName();

const char* GetCircuitChecksum();
CircuitStatistics GetMultiplyCircuitStatistics();
CircuitStatistics GetFFTCircuitStatistics();
CircuitStatistics GetIFFTCircuitStatistics();

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
