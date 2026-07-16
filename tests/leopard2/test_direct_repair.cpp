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
#include <stdint.h>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <vector>

namespace {

using leopard2_test::BinaryField;
using leopard2_test::DirectRepairCrossoverModel;
using leopard2_test::DirectRepairPlan;
using leopard2_test::Element;
using leopard2_test::Matrix;
using leopard2_test::ProfileKind;
using leopard2_test::ProfileLayout;

struct Counts
{
    uint64_t profiles;
    uint64_t plans;
    uint64_t executions;
    uint64_t recovered_symbols;
    uint64_t random_cases;
    uint64_t invalid_cases;

    Counts()
        : profiles(0), plans(0), executions(0), recovered_symbols(0)
        , random_cases(0), invalid_cases(0)
    {
    }
};

class Random
{
public:
    explicit Random(uint64_t seed) : state_(seed) {}

    uint64_t next()
    {
        uint64_t value = state_;
        value ^= value >> 12;
        value ^= value << 25;
        value ^= value >> 27;
        state_ = value;
        return value * UINT64_C(2685821657736338717);
    }

private:
    uint64_t state_;
};

void require(bool condition, const std::string& message)
{
    if (!condition)
        throw std::runtime_error(message);
}

std::vector<std::vector<unsigned> > subsets_of_size(unsigned n, unsigned k)
{
    std::vector<std::vector<unsigned> > result;
    if (k > n)
        return result;
    if (k == 0)
    {
        result.push_back(std::vector<unsigned>());
        return result;
    }
    std::vector<unsigned> subset(k, 0);
    for (unsigned i = 0; i < k; ++i)
        subset[i] = i;
    for (;;)
    {
        result.push_back(subset);
        int i = static_cast<int>(k) - 1;
        while (i >= 0 && subset[static_cast<unsigned>(i)] == n - k + static_cast<unsigned>(i))
            --i;
        if (i < 0)
            break;
        ++subset[static_cast<unsigned>(i)];
        for (unsigned j = static_cast<unsigned>(i) + 1; j < k; ++j)
            subset[j] = subset[j - 1] + 1;
    }
    return result;
}

std::vector<Element> encode(
    const BinaryField& field,
    const Matrix& generator,
    const std::vector<Element>& message)
{
    return leopard2_test::matrix_vector_multiply(field, generator, message);
}

bool plans_equal(const DirectRepairPlan& a, const DirectRepairPlan& b)
{
    if (a.original_count != b.original_count ||
        a.recovery_count != b.recovery_count ||
        a.missing_originals != b.missing_originals ||
        a.selected_parities != b.selected_parities ||
        a.output_terms.size() != b.output_terms.size() ||
        a.setup_counts.inversions != b.setup_counts.inversions ||
        a.setup_counts.multiplications != b.setup_counts.multiplications ||
        a.setup_counts.xors != b.setup_counts.xors ||
        a.execution_inputs_per_symbol != b.execution_inputs_per_symbol ||
        a.execution_multiplications_per_symbol !=
            b.execution_multiplications_per_symbol ||
        a.execution_xors_per_symbol != b.execution_xors_per_symbol)
    {
        return false;
    }
    for (unsigned output = 0; output < a.output_terms.size(); ++output)
    {
        if (a.output_terms[output].size() != b.output_terms[output].size())
            return false;
        for (unsigned i = 0; i < a.output_terms[output].size(); ++i)
        {
            const leopard2_test::DirectRepairTerm& left = a.output_terms[output][i];
            const leopard2_test::DirectRepairTerm& right = b.output_terms[output][i];
            if (left.parity != right.parity || left.index != right.index ||
                left.coefficient != right.coefficient)
                return false;
        }
    }
    return true;
}

void validate_plan_invariants(
    const BinaryField& field,
    const DirectRepairPlan& plan,
    const std::vector<uint8_t>& parity_present)
{
    require(parity_present.size() == plan.recovery_count,
        "plan invariant received a mismatched parity bitmap");
    std::vector<unsigned> expected_parities;
    for (unsigned parity = 0;
         parity < parity_present.size() &&
            expected_parities.size() < plan.missing_originals.size();
         ++parity)
    {
        if (parity_present[parity])
            expected_parities.push_back(parity);
    }
    require(plan.selected_parities == expected_parities,
        "plan did not select the lowest-index available parity rows");
    require(plan.output_terms.size() == plan.missing_originals.size(),
        "plan output count does not match the missing-original list");

    uint64_t input_count = 0;
    uint64_t multiplication_count = 0;
    uint64_t xor_count = 0;
    for (unsigned output = 0; output < plan.output_terms.size(); ++output)
    {
        const std::vector<leopard2_test::DirectRepairTerm>& terms =
            plan.output_terms[output];
        std::vector<uint8_t> seen_original(plan.original_count, 0);
        std::vector<uint8_t> seen_parity(plan.recovery_count, 0);
        input_count += static_cast<uint64_t>(terms.size());
        if (!terms.empty())
            xor_count += static_cast<uint64_t>(terms.size() - 1);
        for (unsigned i = 0; i < terms.size(); ++i)
        {
            const leopard2_test::DirectRepairTerm& term = terms[i];
            require(term.coefficient != 0 && term.coefficient < field.order(),
                "plan retained a zero or out-of-field coefficient");
            if (term.coefficient != 1)
                ++multiplication_count;
            if (term.parity)
            {
                require(term.index < plan.recovery_count &&
                        parity_present[term.index] != 0 &&
                        std::binary_search(plan.selected_parities.begin(),
                            plan.selected_parities.end(), term.index),
                    "plan references an unavailable parity row");
                require(!seen_parity[term.index],
                    "plan contains a duplicate parity term");
                seen_parity[term.index] = 1;
            }
            else
            {
                require(term.index < plan.original_count &&
                        !std::binary_search(plan.missing_originals.begin(),
                            plan.missing_originals.end(), term.index),
                    "plan references a missing original as an input");
                require(!seen_original[term.index],
                    "plan contains a duplicate original term");
                seen_original[term.index] = 1;
            }
        }
    }
    require(input_count == plan.execution_inputs_per_symbol,
        "plan input operation count is incorrect");
    require(multiplication_count == plan.execution_multiplications_per_symbol,
        "plan fixed-multiplication count is incorrect");
    require(xor_count == plan.execution_xors_per_symbol,
        "plan XOR count is incorrect");
}

void validate_execution(
    const BinaryField& field,
    const Matrix& generator,
    const DirectRepairPlan& plan,
    const std::vector<Element>& message,
    Counts* counts)
{
    const std::vector<Element> codeword = encode(field, generator, message);
    std::vector<Element> parities(plan.recovery_count, 0);
    for (unsigned i = 0; i < plan.recovery_count; ++i)
        parities[i] = codeword[plan.original_count + i];
    std::vector<Element> supplied_originals = message;
    for (unsigned i = 0; i < plan.missing_originals.size(); ++i)
        supplied_originals[plan.missing_originals[i]] =
            static_cast<Element>((i + 1) & (field.order() - 1));
    const std::vector<Element> recovered = leopard2_test::execute_direct_repair(
        field, plan, supplied_originals, parities);
    require(recovered.size() == plan.missing_originals.size(),
        "direct repair returned the wrong output count");
    for (unsigned i = 0; i < recovered.size(); ++i)
    {
        if (recovered[i] != message[plan.missing_originals[i]])
        {
            std::ostringstream stream;
            stream << "direct repair mismatch at missing output " << i;
            throw std::runtime_error(stream.str());
        }
        ++counts->recovered_symbols;
    }
    ++counts->executions;
}

void enumerate_messages(
    const BinaryField& field,
    const Matrix& generator,
    const DirectRepairPlan& plan,
    Counts* counts)
{
    const unsigned k = plan.original_count;
    if (k <= 3)
    {
        uint64_t message_count = 1;
        for (unsigned i = 0; i < k; ++i)
            message_count *= field.order();
        for (uint64_t encoded_message = 0; encoded_message < message_count; ++encoded_message)
        {
            uint64_t value = encoded_message;
            std::vector<Element> message(k, 0);
            for (unsigned i = 0; i < k; ++i)
            {
                message[i] = static_cast<Element>(value % field.order());
                value /= field.order();
            }
            validate_execution(field, generator, plan, message, counts);
        }
        return;
    }

    // Direct repair is GF-linear.  Exercising every coordinate-basis value on
    // every systematic input is exhaustive for the linear map without walking
    // all 16^K messages for the larger small-field profiles.
    for (unsigned original = 0; original < k; ++original)
        for (unsigned value = 0; value < field.order(); ++value)
        {
            std::vector<Element> message(k, 0);
            message[original] = static_cast<Element>(value);
            validate_execution(field, generator, plan, message, counts);
        }
}

void test_exhaustive_small_profiles(Counts* counts)
{
    const BinaryField field = leopard2_test::make_gf4();
    const ProfileKind profiles[] = {
        leopard2_test::kLegacyHigh, leopard2_test::kLow
    };
    for (unsigned profile_i = 0; profile_i < 2; ++profile_i)
        for (unsigned k = 1; k <= 5; ++k)
            for (unsigned r = 1; r <= 5; ++r)
            {
                const ProfileLayout layout = leopard2_test::make_profile_layout(
                    profiles[profile_i], k, r);
                if (layout.parent_size > field.order())
                    continue;
                const Matrix generator =
                    leopard2_test::direct_systematic_generator(field, layout);
                ++counts->profiles;
                const unsigned maximum_losses = std::min(k, r);
                for (unsigned losses = 1; losses <= maximum_losses; ++losses)
                {
                    const std::vector<std::vector<unsigned> > missing_sets =
                        subsets_of_size(k, losses);
                    const std::vector<std::vector<unsigned> > parity_sets =
                        subsets_of_size(r, losses);
                    for (unsigned missing_i = 0; missing_i < missing_sets.size(); ++missing_i)
                        for (unsigned parity_i = 0; parity_i < parity_sets.size(); ++parity_i)
                        {
                            std::vector<uint8_t> parity_present(r, 0);
                            for (unsigned i = 0; i < parity_sets[parity_i].size(); ++i)
                                parity_present[parity_sets[parity_i][i]] = 1;
                            const DirectRepairPlan plan =
                                leopard2_test::make_direct_repair_plan(
                                    field, generator, missing_sets[missing_i], parity_present);
                            validate_plan_invariants(field, plan, parity_present);
                            require(plan.selected_parities == parity_sets[parity_i],
                                "direct repair parity selection is not deterministic");
                            enumerate_messages(field, generator, plan, counts);
                            ++counts->plans;
                        }
                }

                // Surplus parity must not change the lowest-index selection.
                std::vector<uint8_t> all_parity(r, 1);
                const std::vector<unsigned> one_missing(1, k - 1);
                const DirectRepairPlan surplus = leopard2_test::make_direct_repair_plan(
                    field, generator, one_missing, all_parity);
                validate_plan_invariants(field, surplus, all_parity);
                require(surplus.selected_parities.size() == 1 &&
                    surplus.selected_parities[0] == 0,
                    "surplus parity did not choose the lowest deterministic row");
            }
}

struct RandomProfile
{
    ProfileKind profile;
    unsigned k;
    unsigned r;
    unsigned iterations;
};

void shuffle_prefix(std::vector<unsigned>* values, unsigned count, Random* random)
{
    for (unsigned i = 0; i < count; ++i)
    {
        const unsigned selected = i + static_cast<unsigned>(
            random->next() % (values->size() - i));
        std::swap((*values)[i], (*values)[selected]);
    }
}

void test_random_field(
    const BinaryField& field,
    const RandomProfile* profiles,
    unsigned profile_count,
    uint64_t seed,
    Counts* counts)
{
    Random random(seed);
    for (unsigned profile_i = 0; profile_i < profile_count; ++profile_i)
    {
        const RandomProfile& configuration = profiles[profile_i];
        const ProfileLayout layout = leopard2_test::make_profile_layout(
            configuration.profile, configuration.k, configuration.r);
        require(layout.parent_size <= field.order(),
            "random direct-repair profile does not fit the selected field");
        const Matrix generator = leopard2_test::direct_systematic_generator(field, layout);
        for (unsigned iteration = 0; iteration < configuration.iterations; ++iteration)
        {
            const unsigned maximum_losses = std::min(
                std::min(configuration.k, configuration.r), 4u);
            const unsigned losses = 1 + static_cast<unsigned>(random.next() % maximum_losses);
            std::vector<unsigned> original_order(configuration.k, 0);
            std::vector<unsigned> parity_order(configuration.r, 0);
            for (unsigned i = 0; i < configuration.k; ++i)
                original_order[i] = i;
            for (unsigned i = 0; i < configuration.r; ++i)
                parity_order[i] = i;
            shuffle_prefix(&original_order, losses, &random);
            shuffle_prefix(&parity_order, losses, &random);
            std::vector<unsigned> missing(
                original_order.begin(), original_order.begin() + losses);
            std::sort(missing.begin(), missing.end());
            std::vector<uint8_t> parity_present(configuration.r, 0);
            for (unsigned i = 0; i < losses; ++i)
                parity_present[parity_order[i]] = 1;
            for (unsigned i = losses; i < configuration.r; ++i)
                if ((random.next() & 3u) == 0)
                    parity_present[parity_order[i]] = 1;

            const DirectRepairPlan plan = leopard2_test::make_direct_repair_plan(
                field, generator, missing, parity_present);
            validate_plan_invariants(field, plan, parity_present);
            require(std::is_sorted(plan.selected_parities.begin(),
                    plan.selected_parities.end()),
                "random plan parity selection is not sorted");
            std::vector<Element> message(configuration.k, 0);
            for (unsigned i = 0; i < configuration.k; ++i)
                message[i] = static_cast<Element>(random.next() & (field.order() - 1));
            validate_execution(field, generator, plan, message, counts);
            ++counts->plans;
            ++counts->random_cases;
        }
    }
}

void test_random_gf8_gf16(Counts* counts)
{
    const RandomProfile gf8_profiles[] = {
        { leopard2_test::kLegacyHigh, 8, 3, 64 },
        { leopard2_test::kLow, 3, 8, 64 },
        { leopard2_test::kLegacyHigh, 31, 7, 32 },
        { leopard2_test::kLow, 7, 31, 32 }
    };
    test_random_field(leopard2_test::make_legacy_gf8(), gf8_profiles,
        sizeof(gf8_profiles) / sizeof(gf8_profiles[0]),
        UINT64_C(0x8f42c6a197), counts);

    const RandomProfile gf16_profiles[] = {
        // Natural GF16 low parent, including the public-count field boundary.
        { leopard2_test::kLow, 2, 256, 48 },
        // Forced-GF16 arithmetic coverage for multi-loss high and low profiles.
        { leopard2_test::kLegacyHigh, 17, 9, 48 },
        { leopard2_test::kLow, 9, 17, 48 },
        // Natural GF16 high parent.  R=1 also exercises the dense K-term plan.
        { leopard2_test::kLegacyHigh, 257, 1, 8 }
    };
    test_random_field(leopard2_test::make_legacy_gf16(), gf16_profiles,
        sizeof(gf16_profiles) / sizeof(gf16_profiles[0]),
        UINT64_C(0x16f00dcafebeef), counts);
}

template<class Function>
void require_throws(Function function, const char* message, Counts* counts)
{
    bool threw = false;
    try
    {
        function();
    }
    catch (const std::exception&)
    {
        threw = true;
    }
    require(threw, message);
    ++counts->invalid_cases;
}

void test_invalid_inputs(Counts* counts)
{
    const BinaryField field = leopard2_test::make_gf4();
    const ProfileLayout layout = leopard2_test::make_profile_layout(
        leopard2_test::kLow, 3, 3);
    const Matrix generator = leopard2_test::direct_systematic_generator(field, layout);
    require_throws([&]() {
        leopard2_test::make_direct_repair_plan(
            field, generator, std::vector<unsigned>(), std::vector<uint8_t>(3, 1));
    }, "empty loss set was accepted", counts);
    require_throws([&]() {
        const std::vector<unsigned> duplicate = { 1, 1 };
        leopard2_test::make_direct_repair_plan(
            field, generator, duplicate, std::vector<uint8_t>(3, 1));
    }, "duplicate loss set was accepted", counts);
    require_throws([&]() {
        const std::vector<unsigned> unsorted = { 2, 0 };
        leopard2_test::make_direct_repair_plan(
            field, generator, unsorted, std::vector<uint8_t>(3, 1));
    }, "unsorted loss set was accepted", counts);
    require_throws([&]() {
        const std::vector<unsigned> missing = { 0, 1 };
        std::vector<uint8_t> one_parity(3, 0);
        one_parity[0] = 1;
        leopard2_test::make_direct_repair_plan(
            field, generator, missing, one_parity);
    }, "insufficient parity was accepted", counts);
    require_throws([&]() {
        leopard2_test::make_direct_repair_plan(
            field, Matrix(), std::vector<unsigned>(1, 0),
            std::vector<uint8_t>(3, 1));
    }, "empty generator was accepted", counts);
    require_throws([&]() {
        Matrix ragged = generator;
        ragged[3].pop_back();
        leopard2_test::make_direct_repair_plan(
            field, ragged, std::vector<unsigned>(1, 0),
            std::vector<uint8_t>(3, 1));
    }, "ragged generator was accepted", counts);
    require_throws([&]() {
        leopard2_test::make_direct_repair_plan(
            field, generator, std::vector<unsigned>(1, 0),
            std::vector<uint8_t>(2, 1));
    }, "mismatched parity bitmap was accepted", counts);
    require_throws([&]() {
        std::vector<uint8_t> invalid_marker(3, 1);
        invalid_marker[1] = 2;
        leopard2_test::make_direct_repair_plan(
            field, generator, std::vector<unsigned>(1, 0), invalid_marker);
    }, "non-Boolean parity marker was accepted", counts);
    require_throws([&]() {
        leopard2_test::make_direct_repair_plan(
            field, generator, std::vector<unsigned>(1, 3),
            std::vector<uint8_t>(3, 1));
    }, "out-of-range missing original was accepted", counts);
    require_throws([&]() {
        Matrix nonsystematic = generator;
        nonsystematic[0][0] = 0;
        leopard2_test::make_direct_repair_plan(
            field, nonsystematic, std::vector<unsigned>(1, 0),
            std::vector<uint8_t>(3, 1));
    }, "nonsystematic generator was accepted", counts);
    require_throws([&]() {
        Matrix invalid_element = generator;
        invalid_element[3][0] = 1;
        invalid_element[3][1] = static_cast<Element>(field.order());
        leopard2_test::make_direct_repair_plan(
            field, invalid_element, std::vector<unsigned>(1, 0),
            std::vector<uint8_t>(3, 1));
    }, "out-of-field generator element was accepted", counts);
    require_throws([&]() {
        Matrix singular = generator;
        singular[3][0] = 0;
        leopard2_test::make_direct_repair_plan(
            field, singular, std::vector<unsigned>(1, 0),
            std::vector<uint8_t>(3, 1));
    }, "singular repair equation was accepted", counts);

    const DirectRepairPlan valid = leopard2_test::make_direct_repair_plan(
        field, generator, std::vector<unsigned>(1, 0),
        std::vector<uint8_t>(3, 1));
    require_throws([&]() {
        leopard2_test::execute_direct_repair(
            field, valid, std::vector<Element>(2, 0),
            std::vector<Element>(3, 0));
    }, "short original execution vector was accepted", counts);
    require_throws([&]() {
        leopard2_test::execute_direct_repair(
            field, valid, std::vector<Element>(3, 0),
            std::vector<Element>(2, 0));
    }, "short parity execution vector was accepted", counts);
    require_throws([&]() {
        ProfileLayout mismatch = layout;
        ++mismatch.recovery_count;
        leopard2_test::model_direct_repair_crossover(valid, mismatch);
    }, "mismatched profile model was accepted", counts);
    require_throws([&]() {
        ProfileLayout irregular = layout;
        irregular.padded_side = 3;
        leopard2_test::model_direct_repair_crossover(valid, irregular);
    }, "non-power-of-two profile model was accepted", counts);
}

void test_specializations_and_immutability(Counts* counts)
{
    const BinaryField field = leopard2_test::make_legacy_gf8();
    const ProfileLayout layout = leopard2_test::make_profile_layout(
        leopard2_test::kLegacyHigh, 5, 1);
    const Matrix generator = leopard2_test::direct_systematic_generator(field, layout);
    const std::vector<unsigned> missing(1, 2);
    const std::vector<uint8_t> parity_present(1, 1);
    const DirectRepairPlan plan = leopard2_test::make_direct_repair_plan(
        field, generator, missing, parity_present);
    validate_plan_invariants(field, plan, parity_present);
    require(plan.output_terms.size() == 1 && plan.output_terms[0].size() == 5,
        "R=1 repair did not retain the parity and four survivors");
    for (unsigned i = 0; i < plan.output_terms[0].size(); ++i)
        require(plan.output_terms[0][i].coefficient == 1,
            "R=1 repair did not specialize a unit coefficient");
    require(plan.execution_inputs_per_symbol == 5 &&
            plan.execution_multiplications_per_symbol == 0 &&
            plan.execution_xors_per_symbol == 4,
        "R=1 specialized operation counts are incorrect");
    require(plan.setup_counts.inversions == 0 &&
            plan.setup_counts.multiplications == 0 &&
            plan.setup_counts.xors == 4,
        "R=1 setup did not specialize zero/unit arithmetic");

    const DirectRepairPlan snapshot = plan;
    const std::vector<Element> message = { 1, 2, 3, 4, 5 };
    validate_execution(field, generator, plan, message, counts);
    validate_execution(field, generator, plan, message, counts);
    require(plans_equal(plan, snapshot),
        "immutable direct repair plan changed during repeated execution");
    ++counts->plans;

    // A nonsingular systematic test matrix is sufficient to exercise zero
    // pruning even though production RS generator matrices are MDS and are
    // generally dense in this representation.
    const BinaryField small_field = leopard2_test::make_gf4();
    Matrix sparse_generator(4, std::vector<Element>(2, 0));
    sparse_generator[0][0] = 1;
    sparse_generator[1][1] = 1;
    sparse_generator[2][0] = 1;
    sparse_generator[2][1] = 1;
    sparse_generator[3][0] = 1;
    const std::vector<unsigned> both_missing = { 0, 1 };
    const std::vector<uint8_t> both_parities(2, 1);
    const DirectRepairPlan sparse = leopard2_test::make_direct_repair_plan(
        small_field, sparse_generator, both_missing, both_parities);
    validate_plan_invariants(small_field, sparse, both_parities);
    require(sparse.output_terms.size() == 2 &&
            sparse.output_terms[0].size() == 1 &&
            sparse.output_terms[1].size() == 2 &&
            sparse.execution_inputs_per_symbol == 3 &&
            sparse.execution_multiplications_per_symbol == 0 &&
            sparse.execution_xors_per_symbol == 1,
        "zero-coefficient pruning operation counts are incorrect");
    const std::vector<Element> sparse_message = { 6, 11 };
    validate_execution(
        small_field, sparse_generator, sparse, sparse_message, counts);
    ++counts->plans;
}

void test_operation_model(Counts* counts)
{
    const BinaryField field = leopard2_test::make_legacy_gf8();
    const ProfileKind kinds[] = {
        leopard2_test::kLegacyHigh, leopard2_test::kLow
    };
    for (unsigned kind_i = 0; kind_i < 2; ++kind_i)
    {
        const unsigned k = kind_i == 0 ? 31 : 7;
        const unsigned r = kind_i == 0 ? 7 : 31;
        const ProfileLayout layout = leopard2_test::make_profile_layout(kinds[kind_i], k, r);
        const Matrix generator = leopard2_test::direct_systematic_generator(field, layout);
        for (unsigned losses = 1; losses <= 4; ++losses)
        {
            std::vector<unsigned> missing(losses, 0);
            for (unsigned i = 0; i < losses; ++i)
                missing[i] = i;
            std::vector<uint8_t> parity_present(r, 1);
            const DirectRepairPlan plan = leopard2_test::make_direct_repair_plan(
                field, generator, missing, parity_present);
            validate_plan_invariants(field, plan, parity_present);
            const DirectRepairCrossoverModel model =
                leopard2_test::model_direct_repair_crossover(plan, layout);
            require(model.direct_operations_per_symbol ==
                plan.execution_multiplications_per_symbol +
                    plan.execution_xors_per_symbol,
                "direct operation model does not match the plan");
            require(model.direct_setup_operations >= plan.setup_counts.inversions,
                "direct setup model dropped inverse operations");
            require(model.specialized_butterflies_per_symbol != 0,
                "transform model unexpectedly has zero work");
            std::cout << "crossover profile="
                      << (kind_i == 0 ? "high" : "low")
                      << " K=" << k << " R=" << r
                      << " L=" << losses
                      << " direct_ops=" << model.direct_operations_per_symbol
                      << " transform_butterflies="
                      << model.specialized_butterflies_per_symbol
                      << " setup_ops=" << model.direct_setup_operations
                      << " symbolic_direct="
                      << (model.direct_execution_below_transform_count ? 1 : 0)
                      << std::endl;
        }
    }
    (void)counts;
}

} // namespace

int main()
{
    try
    {
        Counts counts;
        test_exhaustive_small_profiles(&counts);
        test_random_gf8_gf16(&counts);
        test_invalid_inputs(&counts);
        test_specializations_and_immutability(&counts);
        test_operation_model(&counts);
        std::cout << "direct repair passed: profiles=" << counts.profiles
                  << " plans=" << counts.plans
                  << " executions=" << counts.executions
                  << " recovered_symbols=" << counts.recovered_symbols
                  << " random_cases=" << counts.random_cases
                  << " invalid_cases=" << counts.invalid_cases
                  << std::endl;
        return 0;
    }
    catch (const std::exception& error)
    {
        std::cerr << "direct repair failed: " << error.what() << std::endl;
        return 1;
    }
}
