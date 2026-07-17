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

// Independent algebra checkpoint for Bead leopard-79h.29.  This executable is
// intentionally not linked into the production dispatcher.  It constructs
//
//     Lambda(x) = product_{e in E} (x + omega_e)
//
// with a balanced dense-coefficient product tree, evaluates Lambda away from
// the roots and Lambda' at the roots, and compares those factors with direct
// field algebra, Leopard's active-parent/direct locator implementations, and a
// separately derived R13 Algorithm-8 epsilon interpolation plus shifted-coset
// evaluator.  Neither candidate is linked into production.

#include "LeopardFF8.h"
#include "LeopardFF16.h"

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>

namespace {

static void Fail(const std::string& message)
{
    std::cerr << "locator product-tree experiment failed: " << message
              << std::endl;
    std::exit(1);
}

static void Require(bool condition, const std::string& message)
{
    if (!condition)
        Fail(message);
}

static uint32_t Mix(uint32_t value)
{
    value ^= value >> 16;
    value *= 0x7feb352du;
    value ^= value >> 15;
    value *= 0x846ca68bu;
    return value ^ (value >> 16);
}

struct EpsilonTransformAccounting
{
    uint64_t field_multiplications = 0;
    uint64_t field_additions = 0;
    uint64_t recursive_transform_calls = 0;
    uint64_t transformed_symbol_visits = 0;
    uint64_t fixed_multiplier_uses = 0;
};

struct EpsilonAccounting
{
    uint64_t interpolation_points = 0;
    uint64_t transform_size = 0;
    uint64_t full_parent_special = 0;
    uint64_t direct_prefix_multiplications = 0;
    uint64_t derivative_multiplications = 0;
    uint64_t algorithm8_recursive_calls = 0;
    uint64_t fixed_setup_unique_multipliers = 0;
    uint64_t fixed_setup_unique_subspace_values = 0;
    uint64_t fixed_setup_field_multiplications = 0;
    uint64_t fixed_setup_field_additions = 0;
    uint64_t fixed_setup_inversions = 0;
    uint64_t logical_scratch_elements = 0;
    uint64_t coefficient_elements = 0;
    uint64_t completed_block_elements = 0;
    uint64_t output_writes = 0;
    EpsilonTransformAccounting algorithm8;
    EpsilonTransformAccounting coset_evaluation;
};

struct Accounting
{
    uint64_t product_tree_multiplications = 0;
    uint64_t product_tree_additions = 0;
    uint64_t product_evaluation_multiplications = 0;
    uint64_t product_evaluation_additions = 0;
    uint64_t direct_factor_multiplications = 0;
    uint64_t direct_log_lookups = 0;
    uint64_t direct_modular_additions = 0;
    uint64_t fwht_modular_additions = 0;
    uint64_t fwht_pointwise_multiplications = 0;
    uint64_t peak_tree_elements = 0;
    uint64_t locator_polynomial_elements = 0;
    uint64_t derivative_polynomial_elements = 0;
    EpsilonAccounting epsilon;
};

struct CaseResult
{
    std::string field;
    unsigned field_bits = 0;
    unsigned parent_size = 0;
    unsigned erasure_count = 0;
    uint32_t seed = 0;
    uint64_t compared_entries = 0;
    uint64_t derivative_entries = 0;
    Accounting accounting;
};

class GF4TestField
{
public:
    typedef uint8_t Element;

    GF4TestField()
    {
        for (unsigned i = 0; i < 16; ++i)
        {
            logarithm_[i] = 0;
            exponential_[i] = 0;
        }
        unsigned value = 1;
        for (unsigned exponent = 0; exponent < 15; ++exponent)
        {
            exponential_[exponent] = static_cast<Element>(value);
            logarithm_[value] = exponent;
            value <<= 1;
            if (value & 16)
                value ^= 0x13;
        }
        exponential_[15] = 1;
    }

    const char* Name() const { return "gf2_4"; }
    unsigned Bits() const { return 4; }
    unsigned Order() const { return 16; }
    unsigned Modulus() const { return 15; }
    unsigned ElementBytes() const { return sizeof(Element); }

    Element Multiply(Element a, Element b) const
    {
        if (a == 0 || b == 0)
            return 0;
        return exponential_[(logarithm_[a] + logarithm_[b]) % 15];
    }

    Element Inverse(Element value) const
    {
        Require(value != 0, "GF(2^4) inverse of zero");
        return exponential_[(15 - logarithm_[value]) % 15];
    }

    unsigned Log(Element value) const
    {
        Require(value != 0, "GF(2^4) logarithm of zero");
        return logarithm_[value];
    }

    void ProductionDirect(
        unsigned n,
        const std::vector<uint8_t>& erasures,
        std::vector<unsigned>& output) const
    {
        DirectLogs(*this, n, erasures, output);
    }

    void ProductionWalsh(
        unsigned n,
        const std::vector<uint8_t>& erasures,
        std::vector<unsigned>& output) const
    {
        output.assign(n, 0);
        for (unsigned i = 0; i < n; ++i)
            output[i] = erasures[i] ? 1 : 0;
        Walsh(output, Modulus());

        std::vector<unsigned> kernel(n, 0);
        for (unsigned i = 1; i < n; ++i)
            kernel[i] = Log(static_cast<Element>(i));
        Walsh(kernel, Modulus());
        const unsigned inverse_n = Order() / n;
        for (unsigned i = 0; i < n; ++i)
            output[i] = (output[i] * kernel[i] * inverse_n) % Modulus();
        Walsh(output, Modulus());
    }

private:
    static void Walsh(std::vector<unsigned>& values, unsigned modulus)
    {
        for (size_t distance = 1; distance < values.size(); distance <<= 1)
        {
            for (size_t base = 0; base < values.size(); base += distance << 1)
            {
                for (size_t i = 0; i < distance; ++i)
                {
                    const unsigned a = values[base + i];
                    const unsigned b = values[base + distance + i];
                    values[base + i] = (a + b) % modulus;
                    values[base + distance + i] =
                        (a + modulus - b) % modulus;
                }
            }
        }
    }

    template<class Field>
    static void DirectLogs(
        const Field& field,
        unsigned n,
        const std::vector<uint8_t>& erasures,
        std::vector<unsigned>& output)
    {
        output.assign(n, 0);
        for (unsigned coordinate = 0; coordinate < n; ++coordinate)
        {
            unsigned value = 0;
            for (unsigned erased = 0; erased < n; ++erased)
            {
                if (erasures[erased] && erased != coordinate)
                {
                    value = (value + field.Log(static_cast<Element>(
                        coordinate ^ erased))) % field.Modulus();
                }
            }
            output[coordinate] = value;
        }
    }

    unsigned logarithm_[16];
    Element exponential_[16];
};

struct GF8Field
{
    typedef leopard::ff8::ffe_t Element;
    const char* Name() const { return "gf8"; }
    unsigned Bits() const { return leopard::ff8::kBits; }
    unsigned Order() const { return leopard::ff8::kOrder; }
    unsigned Modulus() const { return leopard::ff8::kModulus; }
    unsigned ElementBytes() const { return sizeof(Element); }
    Element Multiply(Element a, Element b) const
    {
        return leopard::ff8::MultiplyElements(a, b);
    }
    Element Inverse(Element value) const
    {
        return leopard::ff8::InverseElement(value);
    }
    unsigned Log(Element value) const
    {
        return leopard::ff8::ElementLog(value);
    }
    void ProductionDirect(
        unsigned n,
        const std::vector<uint8_t>& erasures,
        std::vector<unsigned>& output) const
    {
        std::vector<Element> temporary(n);
        leopard::ff8::PrepareDecodeDirect(
            n, erasures.data(), temporary.data());
        Convert(temporary, output);
    }
    void ProductionWalsh(
        unsigned n,
        const std::vector<uint8_t>& erasures,
        std::vector<unsigned>& output) const
    {
        std::vector<Element> temporary(n);
        leopard::ff8::PrepareDecodeWalshActive(
            n, erasures.data(), temporary.data());
        Convert(temporary, output);
    }
private:
    static void Convert(
        const std::vector<Element>& source,
        std::vector<unsigned>& destination)
    {
        destination.resize(source.size());
        for (size_t i = 0; i < source.size(); ++i)
            destination[i] = source[i] == leopard::ff8::kModulus ? 0 : source[i];
    }
};

struct GF16Field
{
    typedef leopard::ff16::ffe_t Element;
    const char* Name() const { return "gf16"; }
    unsigned Bits() const { return leopard::ff16::kBits; }
    unsigned Order() const { return leopard::ff16::kOrder; }
    unsigned Modulus() const { return leopard::ff16::kModulus; }
    unsigned ElementBytes() const { return sizeof(Element); }
    Element Multiply(Element a, Element b) const
    {
        return leopard::ff16::MultiplyElements(a, b);
    }
    Element Inverse(Element value) const
    {
        return leopard::ff16::InverseElement(value);
    }
    unsigned Log(Element value) const
    {
        return leopard::ff16::ElementLog(value);
    }
    void ProductionDirect(
        unsigned n,
        const std::vector<uint8_t>& erasures,
        std::vector<unsigned>& output) const
    {
        std::vector<Element> temporary(n);
        leopard::ff16::PrepareDecodeDirect(
            n, erasures.data(), temporary.data());
        Convert(temporary, output);
    }
    void ProductionWalsh(
        unsigned n,
        const std::vector<uint8_t>& erasures,
        std::vector<unsigned>& output) const
    {
        std::vector<Element> temporary(n);
        leopard::ff16::PrepareDecodeWalshActive(
            n, erasures.data(), temporary.data());
        Convert(temporary, output);
    }
private:
    static void Convert(
        const std::vector<Element>& source,
        std::vector<unsigned>& destination)
    {
        destination.resize(source.size());
        for (size_t i = 0; i < source.size(); ++i)
            destination[i] = source[i] == leopard::ff16::kModulus ? 0 : source[i];
    }
};

// Scalar normalized-LCH transforms used only by the epsilon experiment below.
// For V_j = { omega_0, ..., omega_(2^j-1) }, let
//
//     s_j(x) = product_{a in V_j} (x + a).
//
// The additive recurrence
//
//     s_(j+1)(x) = s_j(x) * (s_j(x) + s_j(omega_(2^j)))
//
// derives every skew factor without using Leopard's FFT tables.  Dividing by
// s_j(omega_(2^j)) gives the normalized LCH butterfly multiplier used in R13.
template<class Field>
class ScalarLch
{
public:
    typedef typename Field::Element Element;

    explicit ScalarLch(const Field& field)
        : field_(field)
        , normalizers_(field.Bits())
        , inverse_normalizers_(field.Bits())
    {
        for (unsigned bit = 0; bit < field.Bits(); ++bit)
        {
            Element value = static_cast<Element>(1u << bit);
            for (unsigned lower = 0; lower < bit; ++lower)
            {
                value = field_.Multiply(
                    value, static_cast<Element>(value ^ normalizers_[lower]));
            }
            Require(value != 0, "zero normalized-LCH basis factor");
            normalizers_[bit] = value;
            inverse_normalizers_[bit] = field_.Inverse(value);
        }
    }

    Element Multiplier(
        unsigned dimension,
        unsigned shift,
        EpsilonAccounting& accounting,
        EpsilonTransformAccounting& stage)
    {
        Require(dimension < field_.Bits(), "LCH dimension outside field");
        Require(shift < field_.Order(), "LCH shift outside field");
        ++stage.fixed_multiplier_uses;

        const uint64_t key =
            (static_cast<uint64_t>(dimension) << 32) | shift;
        typename std::map<uint64_t, Element>::const_iterator found =
            multiplier_cache_.find(key);
        if (found != multiplier_cache_.end())
            return found->second;

        if (used_dimensions_.insert(dimension).second)
        {
            // A reusable codec table needs one normalizer and inverse per
            // dimension.  This accounts for the elementary recurrence above,
            // not the already-existing Leopard table construction.
            accounting.fixed_setup_field_multiplications += dimension;
            accounting.fixed_setup_field_additions += dimension;
            ++accounting.fixed_setup_inversions;
        }

        Element numerator = static_cast<Element>(shift);
        for (unsigned lower = 0; lower < dimension; ++lower)
        {
            numerator = field_.Multiply(
                numerator,
                static_cast<Element>(numerator ^ normalizers_[lower]));
            ++accounting.fixed_setup_field_multiplications;
            ++accounting.fixed_setup_field_additions;
        }
        Element result = 0;
        if (numerator != 0)
        {
            result = field_.Multiply(
                numerator, inverse_normalizers_[dimension]);
            ++accounting.fixed_setup_field_multiplications;
        }
        multiplier_cache_[key] = result;
        ++accounting.fixed_setup_unique_multipliers;
        return result;
    }

    Element MultiplyFixed(
        Element multiplier,
        Element value,
        EpsilonTransformAccounting& accounting) const
    {
        if (multiplier == 0 || value == 0)
            return 0;
        if (multiplier == 1)
            return value;
        ++accounting.field_multiplications;
        return field_.Multiply(multiplier, value);
    }

    Element SubspaceValue(
        unsigned dimension,
        unsigned shift,
        EpsilonAccounting& accounting,
        EpsilonTransformAccounting& stage)
    {
        const uint64_t key =
            (static_cast<uint64_t>(dimension) << 32) | shift;
        typename std::map<uint64_t, Element>::const_iterator found =
            subspace_cache_.find(key);
        if (found != subspace_cache_.end())
            return found->second;
        const Element normalized = Multiplier(
            dimension, shift, accounting, stage);
        Element result = 0;
        if (normalized != 0)
        {
            result = field_.Multiply(normalized, normalizers_[dimension]);
            ++accounting.fixed_setup_field_multiplications;
        }
        subspace_cache_[key] = result;
        ++accounting.fixed_setup_unique_subspace_values;
        return result;
    }

private:
    const Field& field_;
    std::vector<Element> normalizers_;
    std::vector<Element> inverse_normalizers_;
    std::map<uint64_t, Element> multiplier_cache_;
    std::map<uint64_t, Element> subspace_cache_;
    std::set<unsigned> used_dimensions_;
};

template<class Field>
static std::vector<typename Field::Element> LchForwardScalar(
    const Field& field,
    ScalarLch<Field>& lch,
    const std::vector<typename Field::Element>& coefficients,
    unsigned shift,
    EpsilonAccounting& total,
    EpsilonTransformAccounting& accounting)
{
    typedef typename Field::Element Element;
    const size_t length = coefficients.size();
    Require(length != 0 && (length & (length - 1)) == 0,
        "non-power-of-two scalar LCH forward");
    Require(length <= field.Order(), "scalar LCH forward exceeds field");
    ++accounting.recursive_transform_calls;
    accounting.transformed_symbol_visits += length;
    if (length == 1)
        return coefficients;

    const size_t half = length / 2;
    unsigned dimension = 0;
    for (size_t value = half; value > 1; value >>= 1)
        ++dimension;
    const Element multiplier = lch.Multiplier(
        dimension, shift, total, accounting);
    std::vector<Element> first(half), second(half);
    for (size_t i = 0; i < half; ++i)
    {
        const Element low = coefficients[i];
        const Element high = coefficients[half + i];
        first[i] = static_cast<Element>(
            low ^ lch.MultiplyFixed(multiplier, high, accounting));
        second[i] = static_cast<Element>(first[i] ^ high);
        accounting.field_additions += 2;
    }
    std::vector<Element> left = LchForwardScalar(
        field, lch, first, shift, total, accounting);
    std::vector<Element> right = LchForwardScalar(
        field, lch, second, shift ^ static_cast<unsigned>(half),
        total, accounting);
    left.insert(left.end(), right.begin(), right.end());
    return left;
}

template<class Field>
static std::vector<typename Field::Element> LchInverseScalar(
    const Field& field,
    ScalarLch<Field>& lch,
    const std::vector<typename Field::Element>& values,
    unsigned shift,
    EpsilonAccounting& total,
    EpsilonTransformAccounting& accounting)
{
    typedef typename Field::Element Element;
    const size_t length = values.size();
    Require(length != 0 && (length & (length - 1)) == 0,
        "non-power-of-two scalar LCH inverse");
    Require(length <= field.Order(), "scalar LCH inverse exceeds field");
    ++accounting.recursive_transform_calls;
    accounting.transformed_symbol_visits += length;
    if (length == 1)
        return values;

    const size_t half = length / 2;
    std::vector<Element> low_values(values.begin(), values.begin() + half);
    std::vector<Element> high_values(values.begin() + half, values.end());
    std::vector<Element> first = LchInverseScalar(
        field, lch, low_values, shift, total, accounting);
    std::vector<Element> second = LchInverseScalar(
        field, lch, high_values, shift ^ static_cast<unsigned>(half),
        total, accounting);

    unsigned dimension = 0;
    for (size_t value = half; value > 1; value >>= 1)
        ++dimension;
    const Element multiplier = lch.Multiplier(
        dimension, shift, total, accounting);
    std::vector<Element> output(length);
    for (size_t i = 0; i < half; ++i)
    {
        const Element upper = static_cast<Element>(first[i] ^ second[i]);
        output[i] = static_cast<Element>(
            first[i] ^ lch.MultiplyFixed(multiplier, upper, accounting));
        output[half + i] = upper;
        accounting.field_additions += 2;
    }
    return output;
}

template<class Element>
struct TangHanResult
{
    std::vector<Element> coefficients;
    std::vector<Element> completed;
};

// Literal scalar form of R13 Appendix-A Algorithm 8.  The first epsilon
// contiguous evaluations determine the unique degree-bounded polynomial;
// coefficients epsilon..q-1 are zero in the normalized LCH basis, and the
// returned q evaluations complete the enclosing power-of-two coset.
template<class Field>
static TangHanResult<typename Field::Element> TangHanAlgorithm8(
    const Field& field,
    ScalarLch<Field>& lch,
    const std::vector<typename Field::Element>& values,
    unsigned count,
    unsigned parent,
    unsigned shift,
    EpsilonAccounting& accounting)
{
    typedef typename Field::Element Element;
    ++accounting.algorithm8_recursive_calls;
    Require(count >= 1 && count <= parent && values.size() == count,
        "invalid Algorithm 8 recursion");
    if (parent == 1)
    {
        TangHanResult<Element> result;
        result.coefficients = values;
        result.completed = values;
        return result;
    }

    const unsigned half = parent / 2;
    if (count <= half)
    {
        TangHanResult<Element> lower = TangHanAlgorithm8(
            field, lch, values, count, half, shift, accounting);
        std::vector<Element> upper = LchForwardScalar(
            field, lch, lower.coefficients, shift ^ half,
            accounting, accounting.algorithm8);
        lower.coefficients.resize(parent, 0);
        lower.completed.insert(
            lower.completed.end(), upper.begin(), upper.end());
        return lower;
    }

    std::vector<Element> first_values(values.begin(), values.begin() + half);
    std::vector<Element> w = LchInverseScalar(
        field, lch, first_values, shift,
        accounting, accounting.algorithm8);
    std::vector<Element> w_prime = LchForwardScalar(
        field, lch, w, shift ^ half,
        accounting, accounting.algorithm8);
    const unsigned tail_count = count - half;
    std::vector<Element> adjusted(tail_count);
    for (unsigned i = 0; i < tail_count; ++i)
    {
        adjusted[i] = static_cast<Element>(values[half + i] ^ w_prime[i]);
        ++accounting.algorithm8.field_additions;
    }
    TangHanResult<Element> upper = TangHanAlgorithm8(
        field, lch, adjusted, tail_count, half, shift ^ half, accounting);

    unsigned dimension = 0;
    for (unsigned value = half; value > 1; value >>= 1)
        ++dimension;
    const Element multiplier = lch.Multiplier(
        dimension, shift, accounting, accounting.algorithm8);
    std::vector<Element> coefficients(parent), completed(parent);
    for (unsigned i = 0; i < half; ++i)
    {
        coefficients[i] = static_cast<Element>(
            w[i] ^ lch.MultiplyFixed(
                multiplier, upper.coefficients[i], accounting.algorithm8));
        coefficients[half + i] = upper.coefficients[i];
        completed[i] = values[i];
        completed[half + i] = static_cast<Element>(
            upper.completed[i] ^ w_prime[i]);
        accounting.algorithm8.field_additions += 2;
    }
    TangHanResult<Element> result;
    result.coefficients.swap(coefficients);
    result.completed.swap(completed);
    return result;
}

template<class Field>
static std::vector<typename Field::Element> MultiplyPolynomials(
    const Field& field,
    const std::vector<typename Field::Element>& left,
    const std::vector<typename Field::Element>& right,
    Accounting& accounting)
{
    typedef typename Field::Element Element;
    std::vector<Element> product(left.size() + right.size() - 1, 0);
    for (size_t i = 0; i < left.size(); ++i)
    {
        for (size_t j = 0; j < right.size(); ++j)
        {
            product[i + j] ^= field.Multiply(left[i], right[j]);
            ++accounting.product_tree_multiplications;
            ++accounting.product_tree_additions;
        }
    }
    return product;
}

template<class Field>
static std::vector<typename Field::Element> BuildLocatorProductTree(
    const Field& field,
    const std::vector<uint8_t>& erasures,
    Accounting& accounting)
{
    typedef typename Field::Element Element;
    std::vector<std::vector<Element> > level;
    for (unsigned coordinate = 0; coordinate < erasures.size(); ++coordinate)
    {
        if (erasures[coordinate])
        {
            std::vector<Element> factor(2);
            factor[0] = static_cast<Element>(coordinate);
            factor[1] = 1;
            level.push_back(factor);
        }
    }
    if (level.empty())
    {
        accounting.peak_tree_elements = 1;
        accounting.locator_polynomial_elements = 1;
        return std::vector<Element>(1, 1);
    }

    uint64_t current_elements = level.size() * 2u;
    accounting.peak_tree_elements = current_elements;
    while (level.size() > 1)
    {
        std::vector<std::vector<Element> > next;
        next.reserve((level.size() + 1) / 2);
        uint64_t next_elements = 0;
        for (size_t i = 0; i < level.size(); i += 2)
        {
            if (i + 1 == level.size())
                next.push_back(level[i]);
            else
                next.push_back(MultiplyPolynomials(
                    field, level[i], level[i + 1], accounting));
            next_elements += next.back().size();
            accounting.peak_tree_elements = std::max(
                accounting.peak_tree_elements,
                current_elements + next_elements);
        }
        level.swap(next);
        current_elements = next_elements;
    }
    accounting.locator_polynomial_elements = level[0].size();
    return level[0];
}

template<class Element>
static std::vector<Element> FormalDerivative(
    const std::vector<Element>& polynomial,
    Accounting& accounting)
{
    if (polynomial.size() <= 1)
    {
        accounting.derivative_polynomial_elements = 1;
        return std::vector<Element>(1, 0);
    }
    std::vector<Element> derivative(polynomial.size() - 1, 0);
    for (size_t degree = 1; degree < polynomial.size(); degree += 2)
        derivative[degree - 1] = polynomial[degree];
    accounting.derivative_polynomial_elements = derivative.size();
    return derivative;
}

template<class Field>
static typename Field::Element EvaluatePolynomial(
    const Field& field,
    const std::vector<typename Field::Element>& polynomial,
    typename Field::Element point,
    Accounting& accounting)
{
    typedef typename Field::Element Element;
    Require(!polynomial.empty(), "empty polynomial");
    Element value = polynomial.back();
    for (size_t i = polynomial.size() - 1; i != 0; --i)
    {
        value = field.Multiply(value, point) ^ polynomial[i - 1];
        ++accounting.product_evaluation_multiplications;
        ++accounting.product_evaluation_additions;
    }
    return value;
}

template<class Field>
static typename Field::Element DirectFactor(
    const Field& field,
    unsigned coordinate,
    const std::vector<uint8_t>& erasures,
    Accounting& accounting)
{
    typedef typename Field::Element Element;
    Element value = 1;
    for (unsigned erased = 0; erased < erasures.size(); ++erased)
    {
        if (!erasures[erased] || erased == coordinate)
            continue;
        value = field.Multiply(
            value, static_cast<Element>(coordinate ^ erased));
        ++accounting.direct_factor_multiplications;
        ++accounting.direct_log_lookups;
        ++accounting.direct_modular_additions;
    }
    return value;
}

static unsigned CeilPowerOfTwo(unsigned value)
{
    Require(value != 0, "zero power-of-two request");
    unsigned result = 1;
    while (result < value)
    {
        Require(result <= (std::numeric_limits<unsigned>::max() >> 1),
            "power-of-two overflow");
        result <<= 1;
    }
    return result;
}

template<class Field>
static std::vector<typename Field::Element> BuildEpsilonLocatorFactors(
    const Field& field,
    unsigned n,
    const std::vector<uint8_t>& erasures,
    const std::vector<typename Field::Element>& expected_locator_values,
    const std::vector<typename Field::Element>& expected_factors,
    Accounting& accounting)
{
    typedef typename Field::Element Element;
    EpsilonAccounting& epsilon = accounting.epsilon;
    const unsigned erasure_count = static_cast<unsigned>(
        std::count(erasures.begin(), erasures.end(), uint8_t(1)));
    Require(erasures.size() == n && expected_locator_values.size() == n &&
            expected_factors.size() == n,
        "epsilon locator input length mismatch");

    if (erasure_count == 0)
    {
        epsilon.interpolation_points = 1;
        epsilon.transform_size = 1;
        epsilon.logical_scratch_elements = 1;
        epsilon.coefficient_elements = 1;
        epsilon.completed_block_elements = 1;
        epsilon.output_writes = n;
        std::vector<Element> result(n, 1);
        Require(result == expected_locator_values && result == expected_factors,
            "empty-erasure locator is not one");
        return result;
    }

    // If every active coordinate is a root, epsilon=E+1 lies outside the
    // active parent (and outside the field when N is the field order).  The
    // active-subspace locator is s_log2(N), whose formal derivative is the
    // same nonzero constant at every root.  Compute that constant once as the
    // direct product at coordinate zero, then scatter it.
    if (erasure_count == n)
    {
        epsilon.interpolation_points = static_cast<uint64_t>(n) + 1;
        epsilon.full_parent_special = 1;
        epsilon.logical_scratch_elements = 1;
        epsilon.output_writes = n;
        Element derivative = 1;
        for (unsigned coordinate = 1; coordinate < n; ++coordinate)
        {
            derivative = field.Multiply(
                derivative, static_cast<Element>(coordinate));
            ++epsilon.derivative_multiplications;
        }
        std::vector<Element> result(n, derivative);
        for (unsigned coordinate = 0; coordinate < n; ++coordinate)
        {
            Require(expected_locator_values[coordinate] == 0,
                "full-parent locator is nonzero at a root");
            Require(result[coordinate] == expected_factors[coordinate],
                "full-parent derivative constant mismatch");
        }
        return result;
    }

    // Normally a monic degree-E locator needs E+1 samples.  When E itself is
    // a power of two q, interpolate only its first q values.  The unique
    // degree-<q interpolant is g(x)=Lambda(x)+s_log2(q)(x), because both
    // monic degree-q polynomials agree on V_q.  On an aligned q-coset the
    // additive subspace polynomial is the constant s_log2(q)(shift), so one
    // XOR correction per output restores Lambda without doubling q.
    const bool power_erasure_count =
        (erasure_count & (erasure_count - 1)) == 0;
    const unsigned interpolation_points = power_erasure_count ?
        erasure_count : erasure_count + 1;
    const unsigned q = CeilPowerOfTwo(erasure_count);
    Require(q <= n && q <= field.Order(),
        "epsilon transform exceeds active parent or field");
    epsilon.interpolation_points = interpolation_points;
    epsilon.transform_size = q;
    epsilon.coefficient_elements = q;
    epsilon.completed_block_elements = q;
    // A reusable in-place realization needs one q-element coefficient block
    // and one q-element value/work block.  The common N-element output is not
    // counted, matching peak_tree_elements above.  The vector-heavy oracle is
    // intentionally not an allocator-RSS claim.
    epsilon.logical_scratch_elements = 2ull * q;

    std::vector<Element> prefix(interpolation_points, 0);
    for (unsigned coordinate = 0;
         coordinate < interpolation_points; ++coordinate)
    {
        if (erasures[coordinate])
            continue;
        Element value = 1;
        for (unsigned erased = 0; erased < n; ++erased)
        {
            if (!erasures[erased])
                continue;
            value = field.Multiply(
                value, static_cast<Element>(coordinate ^ erased));
            ++epsilon.direct_prefix_multiplications;
        }
        prefix[coordinate] = value;
    }

    ScalarLch<Field> lch(field);
    TangHanResult<Element> interpolation = TangHanAlgorithm8(
        field, lch, prefix, interpolation_points, q, 0, epsilon);
    Require(interpolation.coefficients.size() == q &&
            interpolation.completed.size() == q,
        "Algorithm 8 result length mismatch");
    for (unsigned i = interpolation_points; i < q; ++i)
        Require(interpolation.coefficients[i] == 0,
            "Algorithm 8 high coefficient is nonzero");

    // Verify Algorithm 8's two outputs against one another with a fresh LCH
    // context.  This verification work is deliberately excluded from the
    // candidate accounting below.
    EpsilonAccounting verification_accounting;
    ScalarLch<Field> verification_lch(field);
    const std::vector<Element> verified_completion = LchForwardScalar(
        field, verification_lch, interpolation.coefficients, 0,
        verification_accounting, verification_accounting.algorithm8);
    Require(verified_completion == interpolation.completed,
        "Algorithm 8 coefficient/completion mismatch");

    std::vector<Element> raw_locator(n);
    std::copy(interpolation.completed.begin(), interpolation.completed.end(),
        raw_locator.begin());
    epsilon.output_writes += q;
    for (unsigned shift = q; shift < n; shift += q)
    {
        std::vector<Element> block = LchForwardScalar(
            field, lch, interpolation.coefficients, shift,
            epsilon, epsilon.coset_evaluation);
        if (power_erasure_count)
        {
            unsigned dimension = 0;
            for (unsigned value = q; value > 1; value >>= 1)
                ++dimension;
            const Element correction = lch.SubspaceValue(
                dimension, shift, epsilon, epsilon.coset_evaluation);
            if (correction != 0)
            {
                for (unsigned i = 0; i < q; ++i)
                    block[i] ^= correction;
                epsilon.coset_evaluation.field_additions += q;
            }
        }
        std::copy(block.begin(), block.end(), raw_locator.begin() + shift);
        epsilon.output_writes += q;
    }

    std::vector<Element> result = raw_locator;
    for (unsigned coordinate = 0; coordinate < n; ++coordinate)
    {
        if (raw_locator[coordinate] != expected_locator_values[coordinate])
        {
            std::ostringstream message;
            message << field.Name() << " epsilon locator evaluation mismatch n="
                    << n << " erasures=" << erasure_count
                    << " coordinate=" << coordinate;
            Fail(message.str());
        }
        if (!erasures[coordinate])
            continue;

        // This root-only direct product is independent of both the dense
        // formal derivative and the Algorithm 8 interpolation.  The comparison
        // below validates that it equals Lambda'(omega_coordinate).
        Element derivative = 1;
        for (unsigned erased = 0; erased < n; ++erased)
        {
            if (!erasures[erased] || erased == coordinate)
                continue;
            derivative = field.Multiply(
                derivative, static_cast<Element>(coordinate ^ erased));
            ++epsilon.derivative_multiplications;
        }
        result[coordinate] = derivative;
        ++epsilon.output_writes;
    }
    if (result != expected_factors)
    {
        std::ostringstream message;
        message << field.Name() << " epsilon locator factor mismatch n=" << n
                << " erasures=" << erasure_count;
        Fail(message.str());
    }
    return result;
}

static std::vector<uint8_t> MakePattern(
    unsigned n,
    unsigned erasure_count,
    uint32_t seed)
{
    std::vector<unsigned> order(n);
    for (unsigned i = 0; i < n; ++i)
        order[i] = i;
    std::sort(order.begin(), order.end(), [seed](unsigned a, unsigned b) {
        const uint32_t ah = Mix(a ^ seed);
        const uint32_t bh = Mix(b ^ seed);
        return ah == bh ? a < b : ah < bh;
    });
    std::vector<uint8_t> erasures(n, 0);
    for (unsigned i = 0; i < erasure_count; ++i)
        erasures[order[i]] = 1;
    return erasures;
}

template<class Field>
static CaseResult CheckCase(
    const Field& field,
    unsigned n,
    const std::vector<uint8_t>& erasures,
    uint32_t seed)
{
    typedef typename Field::Element Element;
    const unsigned erasure_count = static_cast<unsigned>(
        std::count(erasures.begin(), erasures.end(), uint8_t(1)));
    Accounting accounting;
    const std::vector<Element> locator =
        BuildLocatorProductTree(field, erasures, accounting);
    const std::vector<Element> derivative =
        FormalDerivative(locator, accounting);

    std::vector<unsigned> product_logs(n), independent_direct_logs(n);
    std::vector<Element> locator_values(n), expected_factors(n);
    for (unsigned coordinate = 0; coordinate < n; ++coordinate)
    {
        const Element point = static_cast<Element>(coordinate);
        const Element locator_value =
            EvaluatePolynomial(field, locator, point, accounting);
        if (erasures[coordinate])
            Require(locator_value == 0, "locator root mismatch");

        const Element product_value = erasures[coordinate] ?
            EvaluatePolynomial(field, derivative, point, accounting) :
            locator_value;
        const Element direct_value =
            DirectFactor(field, coordinate, erasures, accounting);
        if (product_value != direct_value)
        {
            std::ostringstream message;
            message << field.Name() << " product/direct mismatch n=" << n
                    << " erasures=" << erasure_count
                    << " coordinate=" << coordinate;
            Fail(message.str());
        }
        Require(product_value != 0, "zero locator factor");
        locator_values[coordinate] = locator_value;
        expected_factors[coordinate] = product_value;
        product_logs[coordinate] = field.Log(product_value);
        independent_direct_logs[coordinate] = field.Log(direct_value);
    }

    std::vector<unsigned> production_direct, production_walsh;
    field.ProductionDirect(n, erasures, production_direct);
    field.ProductionWalsh(n, erasures, production_walsh);
    if (product_logs != independent_direct_logs ||
        product_logs != production_direct ||
        product_logs != production_walsh)
    {
        std::ostringstream message;
        message << field.Name() << " logarithmic differential mismatch n=" << n
                << " erasures=" << erasure_count << " seed=" << seed;
        Fail(message.str());
    }

    const std::vector<Element> epsilon_factors = BuildEpsilonLocatorFactors(
        field, n, erasures, locator_values, expected_factors, accounting);
    Require(epsilon_factors == expected_factors,
        "epsilon locator final differential mismatch");

    unsigned log_n = 0;
    for (unsigned size = n; size > 1; size >>= 1)
        ++log_n;
    accounting.fwht_modular_additions =
        2ull * static_cast<uint64_t>(n) * log_n;
    accounting.fwht_pointwise_multiplications = n;

    CaseResult result;
    result.field = field.Name();
    result.field_bits = field.Bits();
    result.parent_size = n;
    result.erasure_count = erasure_count;
    result.seed = seed;
    result.compared_entries = n * 5ull;
    result.derivative_entries = erasure_count;
    result.accounting = accounting;
    return result;
}

template<class Field>
static void AddDifferentialCases(
    const Field& field,
    bool allow_dense,
    std::vector<CaseResult>& results)
{
    for (unsigned n = 2; n <= field.Order(); n <<= 1)
    {
        std::set<unsigned> counts;
        counts.insert(0);
        counts.insert(1);
        counts.insert(std::min(2u, n));
        counts.insert(std::min(7u, n));
        counts.insert(std::min(14u, n));
        counts.insert(std::min(31u, n));
        counts.insert(std::min(64u, n));
        if (allow_dense || n <= 256)
        {
            counts.insert(n / 4);
            counts.insert(n / 2);
            counts.insert(n);
        }
        unsigned ordinal = 0;
        for (std::set<unsigned>::const_iterator it = counts.begin();
             it != counts.end(); ++it, ++ordinal)
        {
            const uint32_t seed = 0x9e3779b9u ^
                (field.Bits() << 24) ^ (n * 131u) ^ (ordinal * 0x45d9f3bu);
            results.push_back(CheckCase(
                field, n, MakePattern(n, *it, seed), seed));
        }
    }
}

static void AddGF4Exhaustive(std::vector<CaseResult>& results)
{
    const GF4TestField field;
    for (unsigned n = 2; n <= field.Order(); n <<= 1)
    {
        const uint32_t pattern_count = 1u << n;
        for (uint32_t pattern = 0; pattern < pattern_count; ++pattern)
        {
            std::vector<uint8_t> erasures(n);
            for (unsigned i = 0; i < n; ++i)
                erasures[i] = static_cast<uint8_t>((pattern >> i) & 1u);
            results.push_back(CheckCase(field, n, erasures, pattern));
        }
    }
}

static void WriteEpsilonTransformAccounting(
    std::ostream& output,
    const EpsilonTransformAccounting& value)
{
    output << "{\"field_multiplications\":" << value.field_multiplications
           << ",\"field_additions\":" << value.field_additions
           << ",\"recursive_transform_calls\":"
           << value.recursive_transform_calls
           << ",\"transformed_symbol_visits\":"
           << value.transformed_symbol_visits
           << ",\"fixed_multiplier_uses\":"
           << value.fixed_multiplier_uses << "}";
}

static void WriteEpsilonAccounting(
    std::ostream& output,
    const EpsilonAccounting& value)
{
    output << "{\"interpolation_points\":" << value.interpolation_points
           << ",\"transform_size\":" << value.transform_size
           << ",\"full_parent_special\":" << value.full_parent_special
           << ",\"direct_prefix_multiplications\":"
           << value.direct_prefix_multiplications
           << ",\"derivative_multiplications\":"
           << value.derivative_multiplications
           << ",\"algorithm8_recursive_calls\":"
           << value.algorithm8_recursive_calls
           << ",\"fixed_setup_unique_multipliers\":"
           << value.fixed_setup_unique_multipliers
           << ",\"fixed_setup_unique_subspace_values\":"
           << value.fixed_setup_unique_subspace_values
           << ",\"fixed_setup_field_multiplications\":"
           << value.fixed_setup_field_multiplications
           << ",\"fixed_setup_field_additions\":"
           << value.fixed_setup_field_additions
           << ",\"fixed_setup_inversions\":"
           << value.fixed_setup_inversions
           << ",\"logical_scratch_elements\":"
           << value.logical_scratch_elements
           << ",\"coefficient_elements\":" << value.coefficient_elements
           << ",\"completed_block_elements\":"
           << value.completed_block_elements
           << ",\"output_writes\":" << value.output_writes
           << ",\"algorithm8\":";
    WriteEpsilonTransformAccounting(output, value.algorithm8);
    output << ",\"coset_evaluation\":";
    WriteEpsilonTransformAccounting(output, value.coset_evaluation);
    output << "}";
}

static void WriteAccounting(std::ostream& output, const Accounting& value)
{
    output << "{\"product_tree_multiplications\":"
           << value.product_tree_multiplications
           << ",\"product_tree_additions\":" << value.product_tree_additions
           << ",\"product_evaluation_multiplications\":"
           << value.product_evaluation_multiplications
           << ",\"product_evaluation_additions\":"
           << value.product_evaluation_additions
           << ",\"direct_factor_multiplications\":"
           << value.direct_factor_multiplications
           << ",\"direct_log_lookups\":" << value.direct_log_lookups
           << ",\"direct_modular_additions\":"
           << value.direct_modular_additions
           << ",\"fwht_modular_additions\":"
           << value.fwht_modular_additions
           << ",\"fwht_pointwise_multiplications\":"
           << value.fwht_pointwise_multiplications
           << ",\"peak_tree_elements\":" << value.peak_tree_elements
           << ",\"locator_polynomial_elements\":"
           << value.locator_polynomial_elements
           << ",\"derivative_polynomial_elements\":"
           << value.derivative_polynomial_elements
           << ",\"epsilon\":";
    WriteEpsilonAccounting(output, value.epsilon);
    output << "}";
}

static void WriteJson(
    std::ostream& output,
    const std::vector<CaseResult>& results)
{
    uint64_t entries = 0;
    uint64_t derivatives = 0;
    uint64_t gf4_cases = 0;
    uint64_t gf8_cases = 0;
    uint64_t gf16_cases = 0;
    uint64_t accounting_cases = 0;
    uint64_t epsilon_full_parent_special_cases = 0;
    for (size_t i = 0; i < results.size(); ++i)
    {
        entries += results[i].compared_entries;
        derivatives += results[i].derivative_entries;
        epsilon_full_parent_special_cases +=
            results[i].accounting.epsilon.full_parent_special;
        if (results[i].field == "gf2_4")
            ++gf4_cases;
        else if (results[i].field == "gf8")
        {
            ++gf8_cases;
            ++accounting_cases;
        }
        else
        {
            ++gf16_cases;
            ++accounting_cases;
        }
    }

    output << "{\n"
           << "  \"schema\": \"leopard2.locator-construction.algebra.v2\",\n"
           << "  \"timing_evidence\": false,\n"
           << "  \"epsilon_candidate\": \"implemented_nonproduction_correctness_checkpoint_pending_timing\",\n"
           << "  \"summary\": {\"cases\":" << results.size()
           << ",\"gf2_4_exhaustive_cases\":" << gf4_cases
           << ",\"gf8_differential_cases\":" << gf8_cases
           << ",\"gf16_differential_cases\":" << gf16_cases
           << ",\"accounting_cases\":" << accounting_cases
           << ",\"locator_entries_compared\":" << entries
           << ",\"derivative_entries_compared\":" << derivatives
           << ",\"epsilon_cases_cross_checked\":" << results.size()
           << ",\"epsilon_full_parent_special_cases\":"
           << epsilon_full_parent_special_cases
           << "},\n"
           << "  \"gf2_4_exhaustive\": {\"all_subsets\":true,"
              "\"parents\":[{\"n\":2,\"patterns\":4},"
              "{\"n\":4,\"patterns\":16},"
              "{\"n\":8,\"patterns\":256},"
              "{\"n\":16,\"patterns\":65536}]},\n"
           << "  \"cases\": [\n";
    uint64_t emitted = 0;
    for (size_t i = 0; i < results.size(); ++i)
    {
        const CaseResult& result = results[i];
        // The exhaustive GF(2^4) sweep is summarized above.  Repeating all
        // 65,812 individual patterns would inflate this deterministic evidence
        // file by tens of megabytes without adding accounting coverage.
        if (result.field == "gf2_4")
            continue;
        output << "    {\"field\":\"" << result.field
               << "\",\"field_bits\":" << result.field_bits
               << ",\"parent_size\":" << result.parent_size
               << ",\"erasure_count\":" << result.erasure_count
               << ",\"seed\":" << result.seed
               << ",\"compared_entries\":" << result.compared_entries
               << ",\"derivative_entries\":" << result.derivative_entries
               << ",\"element_bytes\":"
               << (result.field_bits <= 8 ? 1 : 2)
               << ",\"accounting\":";
        WriteAccounting(output, result.accounting);
        ++emitted;
        output << "}" << (emitted == accounting_cases ? "\n" : ",\n");
    }
    output << "  ]\n}\n";
}

} // namespace

int main(int argc, char** argv)
{
    Require(argc == 1 || argc == 3,
        "usage: locator_product_tree [--output result.json]");
    Require(argc == 1 || std::string(argv[1]) == "--output",
        "expected --output option");
    Require(leopard::ff8::Initialize(), "GF8 initialization");
    Require(leopard::ff16::Initialize(), "GF16 initialization");

    std::vector<CaseResult> results;
    AddGF4Exhaustive(results);
    AddDifferentialCases(GF8Field(), true, results);
    AddDifferentialCases(GF16Field(), false, results);

    if (argc == 3)
    {
        std::ofstream file(argv[2], std::ios::binary | std::ios::trunc);
        Require(file.good(), "unable to open output file");
        WriteJson(file, results);
        file.flush();
        Require(file.good(), "unable to flush output file");
        file.close();
        Require(!file.fail(), "unable to close output file");
    }
    else
    {
        WriteJson(std::cout, results);
    }

    std::cerr << "locator product-tree algebra passed: cases=" << results.size()
              << std::endl;
    return 0;
}
