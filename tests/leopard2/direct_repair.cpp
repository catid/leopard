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

#include "direct_repair.h"

#include <algorithm>
#include <limits>
#include <stdexcept>

namespace leopard2_test {
namespace {

Element setup_multiply(
    const BinaryField& field,
    Element a,
    Element b,
    DirectRepairSetupCounts* counts)
{
    if (a == 0 || b == 0)
        return 0;
    if (a == 1)
        return b;
    if (b == 1)
        return a;
    ++counts->multiplications;
    return field.multiply(a, b);
}

Matrix invert_counted(
    const BinaryField& field,
    const Matrix& input,
    DirectRepairSetupCounts* counts)
{
    const unsigned n = static_cast<unsigned>(input.size());
    if (n == 0)
        throw std::invalid_argument("direct repair matrix is empty");
    Matrix augmented(n, std::vector<Element>(n * 2, 0));
    for (unsigned row = 0; row < n; ++row)
    {
        if (input[row].size() != n)
            throw std::invalid_argument("direct repair matrix is not square");
        for (unsigned column = 0; column < n; ++column)
            augmented[row][column] = input[row][column];
        augmented[row][n + row] = 1;
    }

    for (unsigned column = 0; column < n; ++column)
    {
        unsigned pivot = column;
        while (pivot < n && augmented[pivot][column] == 0)
            ++pivot;
        if (pivot == n)
            throw std::domain_error("selected direct repair equations are singular");
        if (pivot != column)
            std::swap(augmented[pivot], augmented[column]);

        const Element pivot_value = augmented[column][column];
        if (pivot_value != 1)
        {
            const Element inverse = field.inverse(pivot_value);
            ++counts->inversions;
            for (unsigned j = 0; j < n * 2; ++j)
                augmented[column][j] = setup_multiply(
                    field, augmented[column][j], inverse, counts);
        }

        for (unsigned row = 0; row < n; ++row)
        {
            if (row == column || augmented[row][column] == 0)
                continue;
            const Element factor = augmented[row][column];
            for (unsigned j = 0; j < n * 2; ++j)
            {
                const Element product = setup_multiply(
                    field, factor, augmented[column][j], counts);
                if (product != 0)
                {
                    augmented[row][j] ^= product;
                    ++counts->xors;
                }
            }
        }
    }

    Matrix inverse(n, std::vector<Element>(n, 0));
    for (unsigned row = 0; row < n; ++row)
        for (unsigned column = 0; column < n; ++column)
            inverse[row][column] = augmented[row][n + column];
    return inverse;
}

unsigned integer_log2(unsigned value)
{
    unsigned result = 0;
    while (value > 1)
    {
        if ((value & 1u) != 0)
            throw std::invalid_argument("profile padded side is not a power of two");
        value >>= 1;
        ++result;
    }
    return result;
}

uint64_t checked_add(uint64_t a, uint64_t b)
{
    if (b > std::numeric_limits<uint64_t>::max() - a)
        throw std::overflow_error("direct repair operation count overflow");
    return a + b;
}

uint64_t checked_multiply(uint64_t a, uint64_t b)
{
    if (a != 0 && b > std::numeric_limits<uint64_t>::max() / a)
        throw std::overflow_error("direct repair operation count overflow");
    return a * b;
}

} // namespace

DirectRepairPlan make_direct_repair_plan(
    const BinaryField& field,
    const Matrix& generator,
    const std::vector<unsigned>& missing_originals,
    const std::vector<uint8_t>& parity_present)
{
    if (generator.empty() || generator[0].empty() || missing_originals.empty())
        throw std::invalid_argument("invalid direct repair inputs");
    const unsigned k = static_cast<unsigned>(generator[0].size());
    if (generator.size() < k)
        throw std::invalid_argument("systematic generator has too few rows");
    const unsigned r = static_cast<unsigned>(generator.size()) - k;
    if (parity_present.size() != r || missing_originals.size() > r)
        throw std::invalid_argument("invalid direct repair parity availability");
    for (unsigned row = 0; row < generator.size(); ++row)
    {
        if (generator[row].size() != k)
            throw std::invalid_argument("ragged systematic generator");
        for (unsigned column = 0; column < k; ++column)
            if (generator[row][column] >= field.order())
                throw std::invalid_argument("generator element is outside the field");
    }
    for (unsigned row = 0; row < k; ++row)
        for (unsigned column = 0; column < k; ++column)
            if (generator[row][column] != static_cast<Element>(row == column))
                throw std::invalid_argument("generator is not systematic");

    DirectRepairPlan plan;
    plan.original_count = k;
    plan.recovery_count = r;
    plan.missing_originals = missing_originals;
    std::vector<uint8_t> missing_map(k, 0);
    for (unsigned i = 0; i < missing_originals.size(); ++i)
    {
        const unsigned index = missing_originals[i];
        if (index >= k || missing_map[index])
            throw std::invalid_argument("missing-original list is invalid");
        if (i != 0 && missing_originals[i - 1] >= index)
            throw std::invalid_argument("missing-original list must be sorted");
        missing_map[index] = 1;
    }
    for (unsigned parity = 0;
         parity < r && plan.selected_parities.size() < missing_originals.size();
         ++parity)
    {
        if (parity_present[parity] > 1)
            throw std::invalid_argument("parity-present values must be zero or one");
        if (parity_present[parity])
            plan.selected_parities.push_back(parity);
    }
    if (plan.selected_parities.size() != missing_originals.size())
        throw std::invalid_argument("not enough parity equations for direct repair");
    for (unsigned parity = 0; parity < r; ++parity)
        if (parity_present[parity] > 1)
            throw std::invalid_argument("parity-present values must be zero or one");

    const unsigned losses = static_cast<unsigned>(missing_originals.size());
    Matrix missing_matrix(losses, std::vector<Element>(losses, 0));
    for (unsigned equation = 0; equation < losses; ++equation)
    {
        const std::vector<Element>& generator_row =
            generator[k + plan.selected_parities[equation]];
        for (unsigned missing = 0; missing < losses; ++missing)
            missing_matrix[equation][missing] =
                generator_row[missing_originals[missing]];
    }
    const Matrix inverse = invert_counted(field, missing_matrix, &plan.setup_counts);

    plan.output_terms.resize(losses);
    for (unsigned output = 0; output < losses; ++output)
    {
        std::vector<DirectRepairTerm>& terms = plan.output_terms[output];
        for (unsigned equation = 0; equation < losses; ++equation)
        {
            const Element coefficient = inverse[output][equation];
            if (coefficient != 0)
            {
                DirectRepairTerm term = {
                    true, plan.selected_parities[equation], coefficient
                };
                terms.push_back(term);
            }
        }
        for (unsigned original = 0; original < k; ++original)
        {
            if (missing_map[original])
                continue;
            Element coefficient = 0;
            for (unsigned equation = 0; equation < losses; ++equation)
            {
                const Element product = setup_multiply(
                    field,
                    inverse[output][equation],
                    generator[k + plan.selected_parities[equation]][original],
                    &plan.setup_counts);
                if (product != 0)
                {
                    coefficient ^= product;
                    ++plan.setup_counts.xors;
                }
            }
            if (coefficient != 0)
            {
                DirectRepairTerm term = { false, original, coefficient };
                terms.push_back(term);
            }
        }

        plan.execution_inputs_per_symbol += terms.size();
        if (!terms.empty())
            plan.execution_xors_per_symbol += terms.size() - 1;
        for (unsigned term_i = 0; term_i < terms.size(); ++term_i)
            if (terms[term_i].coefficient != 1)
                ++plan.execution_multiplications_per_symbol;
    }
    return plan;
}

std::vector<Element> execute_direct_repair(
    const BinaryField& field,
    const DirectRepairPlan& plan,
    const std::vector<Element>& originals,
    const std::vector<Element>& parities)
{
    if (originals.size() != plan.original_count ||
        parities.size() != plan.recovery_count ||
        plan.output_terms.size() != plan.missing_originals.size())
    {
        throw std::invalid_argument("direct repair execution vectors have wrong sizes");
    }
    std::vector<Element> result(plan.missing_originals.size(), 0);
    for (unsigned output = 0; output < plan.output_terms.size(); ++output)
    {
        const std::vector<DirectRepairTerm>& terms = plan.output_terms[output];
        bool initialized = false;
        for (unsigned term_i = 0; term_i < terms.size(); ++term_i)
        {
            const DirectRepairTerm& term = terms[term_i];
            const Element input = term.parity
                ? parities[term.index] : originals[term.index];
            const Element product = term.coefficient == 1
                ? input : field.multiply(term.coefficient, input);
            if (!initialized)
            {
                result[output] = product;
                initialized = true;
            }
            else
                result[output] ^= product;
        }
    }
    return result;
}

DirectRepairCrossoverModel model_direct_repair_crossover(
    const DirectRepairPlan& plan,
    const ProfileLayout& layout)
{
    if (layout.original_count != plan.original_count ||
        layout.recovery_count != plan.recovery_count ||
        layout.padded_side == 0 || layout.parent_size == 0)
    {
        throw std::invalid_argument("direct repair plan/profile mismatch");
    }
    DirectRepairCrossoverModel model;
    model.direct_operations_per_symbol = checked_add(
        plan.execution_multiplications_per_symbol,
        plan.execution_xors_per_symbol);
    model.specialized_butterflies_per_symbol = checked_multiply(
        layout.parent_size, integer_log2(layout.padded_side));
    model.direct_setup_operations = checked_add(
        checked_add(plan.setup_counts.inversions, plan.setup_counts.multiplications),
        plan.setup_counts.xors);
    model.direct_execution_below_transform_count =
        model.direct_operations_per_symbol < model.specialized_butterflies_per_symbol;
    return model;
}

} // namespace leopard2_test
