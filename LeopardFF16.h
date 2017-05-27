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


//------------------------------------------------------------------------------
// Fast Walsh-Hadamard Transform (FWHT) (mod kModulus)

// Transform for a variable number of bits (up to kOrder)
void FWHT(ffe_t* data, const unsigned bits);

// Transform specialized for the finite field order
void FWHT(ffe_t data[kOrder]);


//------------------------------------------------------------------------------
// Multiplies

// x[] = y[] * m
void mul_mem_set(
    void * LEO_RESTRICT x, const void * LEO_RESTRICT y,
    ffe_t m, uint64_t bytes);


//------------------------------------------------------------------------------
// FFT Operations

// x[] ^= y[] * m, y[] ^= x[]
void fft_butterfly(
    void * LEO_RESTRICT x, void * LEO_RESTRICT y,
    ffe_t m, uint64_t bytes);

// For i = {0, 1, 2, 3}: x_i[] ^= y_i[] * m, y_i[] ^= x_i[]
void fft_butterfly4(
    void * LEO_RESTRICT x_0, void * LEO_RESTRICT y_0,
    void * LEO_RESTRICT x_1, void * LEO_RESTRICT y_1,
    void * LEO_RESTRICT x_2, void * LEO_RESTRICT y_2,
    void * LEO_RESTRICT x_3, void * LEO_RESTRICT y_3,
    ffe_t m, uint64_t bytes);


//------------------------------------------------------------------------------
// IFFT Operations

// y[] ^= x[], x[] ^= y[] * m
void ifft_butterfly(
    void * LEO_RESTRICT x, void * LEO_RESTRICT y,
    ffe_t m, uint64_t bytes);

// For i = {0, 1, 2, 3}: y_i[] ^= x_i[], x_i[] ^= y_i[] * m
void ifft_butterfly4(
    void * LEO_RESTRICT x_0, void * LEO_RESTRICT y_0,
    void * LEO_RESTRICT x_1, void * LEO_RESTRICT y_1,
    void * LEO_RESTRICT x_2, void * LEO_RESTRICT y_2,
    void * LEO_RESTRICT x_3, void * LEO_RESTRICT y_3,
    ffe_t m, uint64_t bytes);


//------------------------------------------------------------------------------
// Encode

void Encode(
    uint64_t buffer_bytes,
    unsigned original_count,
    unsigned recovery_count,
    unsigned m, // = NextPow2(recovery_count) * 2 = work_count
    void* const * const data,
    void** work); // Size of GetEncodeWorkCount()


//------------------------------------------------------------------------------
// Decode

void Decode(
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


}} // namespace leopard::ff16

#endif // LEO_HAS_FF16
