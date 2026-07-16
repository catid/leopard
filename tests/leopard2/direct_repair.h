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

#include "direct_oracle.h"

#include <stdint.h>

#include <vector>

namespace leopard2_test {

struct DirectRepairTerm
{
    bool parity;
    unsigned index;
    Element coefficient;
};

struct DirectRepairSetupCounts
{
    uint64_t inversions;
    uint64_t multiplications;
    uint64_t xors;

    DirectRepairSetupCounts()
        : inversions(0), multiplications(0), xors(0)
    {
    }
};

// Test-side representation of the immutable portion of a production repair
// plan.  Each output has a fixed list of surviving-original and selected-parity
// terms, so byte-heavy execution does no matrix inversion or coefficient work.
struct DirectRepairPlan
{
    unsigned original_count;
    unsigned recovery_count;
    std::vector<unsigned> missing_originals;
    std::vector<unsigned> selected_parities;
    std::vector<std::vector<DirectRepairTerm> > output_terms;
    DirectRepairSetupCounts setup_counts;
    uint64_t execution_inputs_per_symbol;
    uint64_t execution_multiplications_per_symbol;
    uint64_t execution_xors_per_symbol;

    DirectRepairPlan()
        : original_count(0)
        , recovery_count(0)
        , execution_inputs_per_symbol(0)
        , execution_multiplications_per_symbol(0)
        , execution_xors_per_symbol(0)
    {
    }
};

// Selects the lowest-index available parity equations.  For an MDS systematic
// generator, any L parity rows restricted to L missing systematic columns are
// nonsingular.  Throws if inputs are malformed or fewer than L parities exist.
DirectRepairPlan make_direct_repair_plan(
    const BinaryField& field,
    const Matrix& systematic_generator,
    const std::vector<unsigned>& missing_originals,
    const std::vector<uint8_t>& parity_present);

// Returns recovered values in missing_originals order.  Values at missing
// positions in originals are ignored.  All coefficient work is plan-fixed.
std::vector<Element> execute_direct_repair(
    const BinaryField& field,
    const DirectRepairPlan& plan,
    const std::vector<Element>& originals,
    const std::vector<Element>& parities);

struct DirectRepairCrossoverModel
{
    uint64_t direct_operations_per_symbol;
    uint64_t specialized_butterflies_per_symbol;
    uint64_t direct_setup_operations;
    bool direct_execution_below_transform_count;
};

// Coarse algorithmic model only; it is not a calibrated production threshold.
// The transform side uses N*log2(P or T), matching the specialized decoder's
// stated main-step order.  The direct side counts exact nontrivial fixed
// multiplies and XORs produced by this plan.  Setup and execution stay separate.
DirectRepairCrossoverModel model_direct_repair_crossover(
    const DirectRepairPlan& plan,
    const ProfileLayout& layout);

/*
    Production integration proposal
    -------------------------------

    1. Let codec setup precompute barycentric weights for the profile's parent
       systematic evaluation set.  Plan setup evaluates only L selected parity
       rows against K public originals in O(L*K), rather than constructing the
       test oracle's dense generator.
    2. Before locator/FWHT setup, select the lowest available L parity shards,
       invert their L-by-L missing-column submatrix, and fold every surviving
       systematic contribution into the immutable term lists represented here.
    3. Store coefficients in backend-ready fixed-multiplier form.  Execution is
       L output accumulations over selected parity and surviving original input
       shards, with coefficient 0 removed and coefficient 1 specialized to XOR.
       It requires no allocation, locator work, or N-shard transform scratch.
    4. Dispatch using calibrated cells over L, K, shard bytes, field, backend,
       batch size, and reuse count.  This symbolic model only identifies cells
       worth benchmarking; it must not become the production threshold itself.
    5. Keep output alias rules identical to normal decode and build each missing
       output in its own caller/scratch buffer before scatter.  Re-encoding
       parity remains an explicit separate operation.
*/

} // namespace leopard2_test
