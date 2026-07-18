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

#include <stdint.h>

namespace leopard { namespace ff8xor { namespace avx2 {


// These entry points live in a translation unit built with AVX2 enabled.  The
// baseline dispatcher must check both KernelsBuilt() and CPU/OS support before
// entering any processing function.
bool KernelsBuilt();

// Plane operations return the first unprocessed byte offset in each plane.
// start_offset lets an AVX-512 ZMM kernel hand its 32-byte tail to AVX2 without
// repeating work.  The coefficient/skew dispatch has already happened in the
// caller, so these templates add no dispatch inside the payload loop.
template <unsigned Coefficient>
uint64_t Multiply(
    void* destination,
    const void* source,
    uint64_t buffer_bytes,
    uint64_t start_offset);

template <unsigned Skew, bool Inverse>
uint64_t Butterfly(
    void* x,
    void* y,
    uint64_t buffer_bytes,
    uint64_t start_offset);

// Contiguous XOR is used by butterfly skew 255.  Its offset addresses the
// entire buffer instead of one plane.
uint64_t Xor(
    void* destination,
    const void* source,
    uint64_t buffer_bytes,
    uint64_t start_offset);

// Batched contiguous XOR keeps independent streams in one vector loop so the
// core can overlap their load/XOR/store chains.  Pointers are named explicitly
// to prevent a dynamically indexed hot-path register array.
uint64_t Xor2(
    void* destination0,
    const void* source0,
    void* destination1,
    const void* source1,
    uint64_t buffer_bytes,
    uint64_t start_offset);

uint64_t Xor4(
    void* destination0,
    const void* source0,
    void* destination1,
    const void* source1,
    void* destination2,
    const void* source2,
    void* destination3,
    const void* source3,
    uint64_t buffer_bytes,
    uint64_t start_offset);

// One source-arity-specialized row of the combined formal-derivative/top-FFT
// boundary.  The arity switch occurs before the YMM payload loop.
uint64_t FormalDerivativeBoundaryRow(
    unsigned extra_count,
    void* left,
    void* right,
    const void* extra0,
    const void* extra1,
    const void* extra2,
    const void* extra3,
    const void* extra4,
    const void* extra5,
    const void* extra6,
    uint64_t buffer_bytes,
    uint64_t start_offset);


}}} // namespace leopard::ff8xor::avx2
