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

/*
    8-bit Finite Field Math

    This finite field contains 256 elements and so each element is one byte.
    This library is designed for data that is a multiple of 64 bytes in size.
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


//------------------------------------------------------------------------------
// Fast Walsh-Hadamard Transform (FWHT) (mod kModulus)

// Define this to enable the optimized version of FWHT()
#define LEO_FF8_FWHT_OPTIMIZED

// Transform for a variable number of bits (up to kOrder)
void FWHT(ffe_t* data, const unsigned bits);

// Transform specialized for the finite field order
void FWHT(ffe_t data[kOrder]);


//------------------------------------------------------------------------------
// XOR Memory

// x[] ^= y[]
void xor_mem(
    void * LEO_RESTRICT x, const void * LEO_RESTRICT y,
    unsigned bytes);

// For i = {0, 1}: x_i[] ^= x_i[]
void xor_mem2(
    void * LEO_RESTRICT x_0, const void * LEO_RESTRICT y_0,
    void * LEO_RESTRICT x_1, const void * LEO_RESTRICT y_1,
    unsigned bytes);

// For i = {0, 1, 2}: x_i[] ^= x_i[]
void xor_mem3(
    void * LEO_RESTRICT x_0, const void * LEO_RESTRICT y_0,
    void * LEO_RESTRICT x_1, const void * LEO_RESTRICT y_1,
    void * LEO_RESTRICT x_2, const void * LEO_RESTRICT y_2,
    unsigned bytes);


//------------------------------------------------------------------------------
// Multiplies

// x[] = y[] * m
void mul_mem_set(
    void * LEO_RESTRICT x, const void * LEO_RESTRICT y,
    ffe_t m, unsigned bytes);

// For i = {0, 1}: x_i[] *= m
void mul_mem2_inplace(
    void * LEO_RESTRICT x_0,
    void * LEO_RESTRICT x_1,
    ffe_t m, unsigned bytes);


//------------------------------------------------------------------------------
// FFT Operations

// x[] ^= y[] * m, y[] ^= x[]
void mul_fft(
    void * LEO_RESTRICT x, void * LEO_RESTRICT y,
    ffe_t m, unsigned bytes);

// For i = {0, 1}: x_i[] ^= y_i[] * m, y_i[] ^= x_i[]
void mul_fft2(
    void * LEO_RESTRICT x_0, void * LEO_RESTRICT y_0,
    void * LEO_RESTRICT x_1, void * LEO_RESTRICT y_1,
    ffe_t m, unsigned bytes);

// For i = {0, 1, 2}: x_i[] ^= y_i[] * m, y_i[] ^= x_i[]
void mul_fft3(
    void * LEO_RESTRICT x_0, void * LEO_RESTRICT y_0,
    void * LEO_RESTRICT x_1, void * LEO_RESTRICT y_1,
    void * LEO_RESTRICT x_2, void * LEO_RESTRICT y_2,
    ffe_t m, unsigned bytes);


//------------------------------------------------------------------------------
// IFFT Operations

// y[] ^= x[], x[] ^= y[] * m
void mul_ifft(
    void * LEO_RESTRICT x, void * LEO_RESTRICT y,
    ffe_t m, unsigned bytes);

// For i = {0, 1}: y_i[] ^= x_i[], x_i[] ^= y_i[] * m
void mul_ifft2(
    void * LEO_RESTRICT x_0, void * LEO_RESTRICT y_0,
    void * LEO_RESTRICT x_1, void * LEO_RESTRICT y_1,
    ffe_t m, unsigned bytes);

// For i = {0, 1, 2}: y_i[] ^= x_i[], x_i[] ^= y_i[] * m
void mul_ifft3(
    void * LEO_RESTRICT x_0, void * LEO_RESTRICT y_0,
    void * LEO_RESTRICT x_1, void * LEO_RESTRICT y_1,
    void * LEO_RESTRICT x_2, void * LEO_RESTRICT y_2,
    ffe_t m, unsigned bytes);


//------------------------------------------------------------------------------
// API

// Returns false if the self-test fails
bool Initialize();


}} // namespace leopard::ff8
