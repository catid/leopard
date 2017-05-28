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
// Fast Walsh-Hadamard Transform (FWHT) (mod kModulus)

// Transform for a variable number of bits (up to kOrder)
void FWHT(ffe_t* data, const unsigned bits);

// Transform specialized for the finite field order
void FWHT(ffe_t data[kOrder]);


//------------------------------------------------------------------------------
// Multiplies

// x[] = exp(log(y[]) + log_m)
void mul_mem(
    void * LEO_RESTRICT x, const void * LEO_RESTRICT y,
    ffe_t log_m, uint64_t bytes);


//------------------------------------------------------------------------------
// FFT Operations

/*
    Precondition: log_m != kModulus

    x[] ^= exp(log(y[]) + log_m)
    y[] ^= x[]
*/
void fft_butterfly(
    void * LEO_RESTRICT x, void * LEO_RESTRICT y,
    ffe_t log_m, uint64_t bytes);

#ifdef LEO_USE_VECTOR4_OPT

// Unroll 4 rows at a time
void fft_butterfly4(
    void * LEO_RESTRICT x_0, void * LEO_RESTRICT y_0,
    void * LEO_RESTRICT x_1, void * LEO_RESTRICT y_1,
    void * LEO_RESTRICT x_2, void * LEO_RESTRICT y_2,
    void * LEO_RESTRICT x_3, void * LEO_RESTRICT y_3,
    ffe_t log_m, uint64_t bytes);

#endif // LEO_USE_VECTOR4_OPT


//------------------------------------------------------------------------------
// IFFT Operations

/*
    Precondition: log_m != kModulus

    y[] ^= x[]
    x[] ^= exp(log(y[]) + log_m)
*/
void ifft_butterfly(
    void * LEO_RESTRICT x, void * LEO_RESTRICT y,
    ffe_t log_m, uint64_t bytes);

#ifdef LEO_USE_VECTOR4_OPT

// Unroll 4 rows at a time
void ifft_butterfly4(
    void * LEO_RESTRICT x_0, void * LEO_RESTRICT y_0,
    void * LEO_RESTRICT x_1, void * LEO_RESTRICT y_1,
    void * LEO_RESTRICT x_2, void * LEO_RESTRICT y_2,
    void * LEO_RESTRICT x_3, void * LEO_RESTRICT y_3,
    ffe_t log_m, uint64_t bytes);

#endif // LEO_USE_VECTOR4_OPT


//------------------------------------------------------------------------------
// FFT

/*
    if (log_m != kModulus)
        x[] ^= exp(log(y[]) + log_m)
    y[] ^= x[]
*/
void VectorFFTButterfly(
    const uint64_t bytes,
    unsigned count,
    void** x,
    void** y,
    const ffe_t log_m);

/*
    y[] ^= x[]
    if (log_m != kModulus)
        x[] ^= exp(log(y[]) + log_m)
*/
void VectorIFFTButterfly(
    const uint64_t bytes,
    unsigned count,
    void** x,
    void** y,
    const ffe_t log_m);


//------------------------------------------------------------------------------
// Reed-Solomon Encode

void ReedSolomonEncode(
    uint64_t buffer_bytes,
    unsigned original_count,
    unsigned recovery_count,
    unsigned m, // = NextPow2(recovery_count) * 2 = work_count
    void* const * const data,
    void** work); // Size of GetEncodeWorkCount()


//------------------------------------------------------------------------------
// Reed-Solomon Decode

void ReedSolomonDecode(
    uint64_t buffer_bytes,
    unsigned original_count,
    unsigned recovery_count,
    unsigned m, // = NextPow2(recovery_count)
    unsigned n, // = NextPow2(m + original_count) = work_count
    void* const * const original, // original_count entries
    void* const * const recovery, // recovery_count entries
    void** work); // n entries


//------------------------------------------------------------------------------
// API

// Returns false if the self-test fails
bool Initialize();


}} // namespace leopard::ff8

#endif // LEO_HAS_FF8
