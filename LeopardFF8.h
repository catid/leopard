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

/*
    8-bit Finite Field Math

    This finite field contains 256 elements and so each element is one byte.
    This library is designed for data that is a multiple of 64 bytes in size.

    Algorithms are described in LeopardCommon.h
*/

namespace leopard { namespace ff8 {


//------------------------------------------------------------------------------
// Datatypes and Constants

// Finite field element type
typedef uint8_t ffe_t;

// Number of bits per element
static const unsigned kBits = 8;

// Finite field order: Number of elements in the field
static const unsigned kOrder = 256;

// Modulus for field operations
static const ffe_t kModulus = 255;

// LFSR Polynomial that generates the field elements
static const unsigned kPolynomial = 0x11D;


//------------------------------------------------------------------------------
// API

// Returns false if the self-test fails
bool Initialize();

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


}} // namespace leopard::ff8

#endif // LEO_HAS_FF8
