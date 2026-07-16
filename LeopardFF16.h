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

#ifdef LEO_HAS_FF16

/*
    16-bit Finite Field Math

    This finite field contains 65536 elements and so each element is one byte.
    This library is designed for data that is a multiple of 64 bytes in size.

    Algorithms are described in LeopardCommon.h
*/

namespace leopard { namespace ff16 {


//------------------------------------------------------------------------------
// Datatypes and Constants

// Finite field element type
typedef uint16_t ffe_t;

// Number of bits per element
static const unsigned kBits = 16;

// Finite field order: Number of elements in the field
static const unsigned kOrder = 65536;

// Modulus for field operations
static const ffe_t kModulus = 65535;

// LFSR Polynomial that generates the field elements
static const unsigned kPolynomial = 0x1002D;


//------------------------------------------------------------------------------
// API

// Returns false if the self-test fails
bool Initialize();

// Internal scalar field helpers used while constructing immutable Leopard2
// direct-repair plans.  Elements use Leopard's legacy Cantor coordinates.
ffe_t MultiplyElements(ffe_t a, ffe_t b);
ffe_t InverseElement(ffe_t value);
ffe_t ElementLog(ffe_t value); // value must be nonzero

// Fixed-multiplier execution helpers.  multiplier_log is produced by
// ElementLog(), source and destination must not overlap, and byte_count must
// contain complete GF16 symbols.  Full 64-byte ALTMAP tiles use the active SIMD
// backend.  A final compact tile stores q low bytes followed by q high bytes.
void MultiplyBytes(
    void* destination,
    const void* source,
    ffe_t multiplier_log,
    uint64_t byte_count);
void MultiplyAddBytes(
    void* destination,
    const void* source,
    ffe_t multiplier_log,
    uint64_t byte_count);

void ReedSolomonEncode(
    uint64_t buffer_bytes,
    unsigned original_count,
    unsigned recovery_count,
    unsigned m, // = NextPow2(recovery_count)
    const void* const * const data,
    void** work); // m * 2 elements

// Encodes the low-rate coordinate profile:
//   [ data: 0 .. p-1 ][ recovery: p .. ]
// original_count may be smaller than p; the remaining data coordinates are
// shortened to zero.  Null recovery output pointers are not written.  p must
// be a power of two in [2, kOrder / 2].
void ReedSolomonEncodeLow(
    uint64_t buffer_bytes,
    unsigned original_count,
    unsigned recovery_count,
    unsigned p, // = NextPow2(original_count)
    const void* const * const data,
    void* const * const recovery,
    void** work); // p * 2 elements

// Precomputes logarithmic error-locator evaluations for an active power-of-two
// parent.  erasures contains n bytes, one for each active coordinate.
void PrepareDecode(
    unsigned n,
    const uint8_t* erasures,
    ffe_t* locator_logs); // n elements

// Reuses locator contributions for permanent erasures (for example punctured
// parent coordinates) and adds only the pattern-specific erasures.  The base
// array must have been produced by PrepareDecode for permanent_erasures.
void PrepareDecodeWithPermanent(
    unsigned n,
    const uint8_t* erasures,
    const uint8_t* permanent_erasures,
    const ffe_t* permanent_locator_logs,
    ffe_t* locator_logs); // n elements

// Full-field FWHT implementation retained as an independent differential
// oracle for the sparse active-parent locator path.
void PrepareDecodeWalshReference(
    unsigned n,
    const uint8_t* erasures,
    ffe_t* locator_logs); // n elements

// Returns the deterministic built-in locator setup decision for this field.
bool IsDirectLocatorPreferred(unsigned n, unsigned erasure_count);

// Generic active-parent decoder using locator evaluations from PrepareDecode.
// Null coordinate_data pointers are treated as zero.  Only coordinates marked
// in requested_outputs are revealed, in-place, in the corresponding work slot.
void ReedSolomonDecodePrepared(
    uint64_t buffer_bytes,
    unsigned n,
    const void* const * const coordinate_data, // n elements
    const uint8_t* requested_outputs, // n elements
    const ffe_t* locator_logs, // n elements
    void** work); // n elements

// Precompute the active-parent fixed factors used by the IT2026 specialized
// low/high original-recovery decoders.
void PrepareLowDecode(
    unsigned n,
    unsigned p,
    ffe_t* block_factors); // n / p - 1 elements

void PrepareHighDecode(
    unsigned n,
    unsigned t,
    ffe_t* output_factors); // n elements; t..n-1 are initialized

#if defined(LEO2_ENABLE_TEST_HOOKS)
// Internal test-only access to the exact production LCH kernels and constants.
// These are deliberately absent from leopard2.h and from non-test builds.
ffe_t TestOnlyFFTMultiplier(unsigned logical_index);
ffe_t TestOnlySubspaceDerivative(unsigned size);
ffe_t TestOnlySubspaceAt(unsigned size, unsigned shift);
ffe_t TestOnlyLchNormalizer(unsigned index);
void TestOnlyLchForward(
    uint64_t buffer_bytes,
    unsigned size,
    unsigned shift,
    unsigned requested_output_count,
    void** work);
void TestOnlyLchInverse(
    uint64_t buffer_bytes,
    unsigned size,
    unsigned shift,
    unsigned known_input_count,
    void** work);
void TestOnlyAddFormalDerivative(
    uint64_t buffer_bytes,
    unsigned size,
    void** work);
#endif

void ReedSolomonDecodeLowPrepared(
    uint64_t buffer_bytes,
    unsigned n,
    unsigned p,
    const void* const * const coordinate_data,
    const uint8_t* requested_outputs,
    const ffe_t* locator_logs,
    const ffe_t* block_factors,
    void** work);

void ReedSolomonDecodeHighPrepared(
    uint64_t buffer_bytes,
    unsigned n,
    unsigned t,
    const void* const * const coordinate_data,
    const uint8_t* requested_outputs,
    const ffe_t* locator_logs,
    const ffe_t* output_factors,
    void** work);

void ReedSolomonDecode(
    uint64_t buffer_bytes,
    unsigned original_count,
    unsigned recovery_count,
    unsigned m, // = NextPow2(recovery_count)
    unsigned n, // = NextPow2(m + original_count)
    const void* const * const original, // original_count elements
    const void* const * const recovery, // recovery_count elements
    void** work); // n elements


}} // namespace leopard::ff16

#endif // LEO_HAS_FF16
