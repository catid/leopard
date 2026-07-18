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

namespace leopard { namespace ff8xor { namespace avx512 {


// These entry points live in a translation unit built with AVX-512 enabled.
// The baseline dispatcher must check both KernelsBuilt() and the CPU/OS state
// before calling any of the processing functions.
bool KernelsBuilt();

// Plane operations return the number of bytes processed in each plane.  The
// caller completes the AVX2, 128-bit, and uint64_t tails starting at that
// offset.  The caller's existing coefficient/skew dispatch selects one of
// these direct specializations; they add no second runtime dispatch.
template <unsigned Coefficient>
uint64_t Multiply256(
    void* destination,
    const void* source,
    uint64_t buffer_bytes);

template <unsigned Coefficient>
uint64_t Multiply512(
    void* destination,
    const void* source,
    uint64_t buffer_bytes);

template <unsigned Skew, bool Inverse>
uint64_t Butterfly256(
    void* x,
    void* y,
    uint64_t buffer_bytes);

template <unsigned Skew, bool Inverse>
uint64_t Butterfly512(
    void* x,
    void* y,
    uint64_t buffer_bytes);

// Contiguous XOR is used by butterfly skew 255.  These return the total byte
// offset rather than a per-plane offset because the sentinel operation does
// not need to address the planes independently.
uint64_t Xor256(
    void* destination,
    const void* source,
    uint64_t buffer_bytes);

uint64_t Xor512(
    void* destination,
    const void* source,
    uint64_t buffer_bytes);

uint64_t Xor2_256(
    void* destination0,
    const void* source0,
    void* destination1,
    const void* source1,
    uint64_t buffer_bytes);

uint64_t Xor2_512(
    void* destination0,
    const void* source0,
    void* destination1,
    const void* source1,
    uint64_t buffer_bytes);

uint64_t Xor4_256(
    void* destination0,
    const void* source0,
    void* destination1,
    const void* source1,
    void* destination2,
    const void* source2,
    void* destination3,
    const void* source3,
    uint64_t buffer_bytes);

uint64_t Xor4_512(
    void* destination0,
    const void* source0,
    void* destination1,
    const void* source1,
    void* destination2,
    const void* source2,
    void* destination3,
    const void* source3,
    uint64_t buffer_bytes);


}}} // namespace leopard::ff8xor::avx512
