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

namespace leopard { namespace ff8xor { namespace transpose {


// These boundary helpers convert between ordinary packed FF8 shards and the
// experimental native eight-plane representation.  They allocate no memory.
// Source and destination must not overlap, and buffer_bytes must be divisible
// by eight.  Native codec payload processing never calls these helpers.
enum class Mode
{
    Auto,
    Portable,
    Avx2,
    Avx512Bitalg,
    Avx512Vbmi
};

bool IsModeAvailable(Mode mode);

// Pure predicates used by runtime dispatch and exhaustive synthetic-feature
// tests.  The BITALG and VBMI objects intentionally have different feature
// contracts and neither requires AVX-512VL.
bool IsAVX512BitalgSupported(
    uint32_t maximum_basic_leaf,
    uint32_t leaf1_ecx,
    uint32_t leaf7_ebx,
    uint32_t leaf7_ecx,
    uint64_t xcr0);

bool IsAVX512VbmiSupported(
    uint32_t maximum_basic_leaf,
    uint32_t leaf1_ecx,
    uint32_t leaf7_ebx,
    uint32_t leaf7_ecx,
    uint64_t xcr0);

void PackedToPlane(
    const void* packed,
    void* plane,
    uint64_t buffer_bytes,
    Mode mode = Mode::Auto);

// The same blocked dispatch is available in the inverse direction.  A
// portable word-transpose tail/reference remains available for every target.
void PlaneToPacked(
    const void* plane,
    void* packed,
    uint64_t buffer_bytes,
    Mode mode = Mode::Auto);


}}} // namespace leopard::ff8xor::transpose
