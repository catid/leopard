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

#include "leopard2.h"

#include <stddef.h>
#include <stdint.h>

namespace leopard2_internal {

/*
    Pure policy predicate kept separate from execution so every boundary of the
    evidence-scoped automatic crossover can be tested without timing inference.
*/
static inline bool ShouldUseBalancedGenericDecode(
    leo2_profile profile,
    leo2_field field,
    uint32_t original_count,
    uint32_t recovery_count,
    uint32_t padded_side,
    uint32_t parent_count,
    uint32_t missing_original_count,
    size_t rounded_shard_bytes,
    leo2_backend backend)
{
    return profile == LEO2_PROFILE_LEGACY_HIGH_V1 &&
        field == LEO2_FIELD_GF8 &&
        original_count == 128 && recovery_count == 128 &&
        padded_side == 128 && parent_count == 256 &&
        missing_original_count == 128 &&
        rounded_shard_bytes >= 256 && rounded_shard_bytes <= 1024 * 1024 &&
        (backend == LEO2_BACKEND_SCALAR ||
         backend == LEO2_BACKEND_SSSE3 ||
         backend == LEO2_BACKEND_AVX2);
}

/*
    Offline-calibrated exception to the side-sized high decoder.  Both kernels
    implement the same legacy-high mathematics and wire profile; this policy
    selects the retained regular N-slot traversal only where three-round pinned
    evidence found a credible tiled regression.  Batch calls are handled by the
    caller because two or more stripes restore the tiled cache advantage.
*/
static inline bool ShouldUseMaterializedHighDecode(
    leo2_profile profile,
    leo2_field field,
    uint32_t original_count,
    uint32_t recovery_count,
    uint32_t padded_side,
    uint32_t parent_count,
    uint32_t missing_original_count,
    size_t rounded_shard_bytes,
    leo2_backend backend)
{
    if (profile != LEO2_PROFILE_LEGACY_HIGH_V1 ||
        field != LEO2_FIELD_GF8 ||
        original_count != 224 || recovery_count != 32 ||
        padded_side != 32 || parent_count != 256 ||
        missing_original_count == 0 || missing_original_count > 8 ||
        rounded_shard_bytes > 64 * 1024)
        return false;

    if (backend == LEO2_BACKEND_AVX2)
        return rounded_shard_bytes >= 24 * 1024;
    if (backend == LEO2_BACKEND_SSSE3)
        return rounded_shard_bytes >= 32 * 1024;
    return false;
}

} // namespace leopard2_internal
