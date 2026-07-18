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

namespace leopard { namespace ff8xor { namespace avx512four {


// This object is compiled separately with AVX-512 enabled.  The baseline
// dispatcher must prove the CPU/OS feature contract before calling Apply512.
bool KernelsBuilt();

unsigned GetTupleCount();
bool GetTuple(
    unsigned tuple_index,
    uint8_t& skew01,
    uint8_t& skew23,
    uint8_t& skew02);

// Returns the generated specialization for a reachable, complete two-layer
// radix-4 unit, or -1 when the tuple is not in the checked corpus.
int FindTupleIndex(uint8_t skew01, uint8_t skew23, uint8_t skew02);

// Apply one complete two-layer transform to four native plane-layout buffers.
// The coefficient/function-table dispatch happens once here, before the
// generated specialization enters its ZMM offset loop.  False means that the
// caller must preserve the two-way schedule (unsupported tuple or narrow tail).
bool Apply512(
    unsigned tuple_index,
    void* a,
    void* b,
    void* c,
    void* d,
    uint64_t buffer_bytes,
    bool inverse,
    bool use_xor3);


}}} // namespace leopard::ff8xor::avx512four
