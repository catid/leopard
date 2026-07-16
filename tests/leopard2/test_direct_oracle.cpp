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

#include "direct_oracle.h"

#include <iostream>
#include <sstream>
#include <stdexcept>

namespace {

using namespace leopard2_test;

struct Counts
{
    uint64_t field_checks;
    uint64_t subspace_checks;
    uint64_t normalization_checks;
    uint64_t high_basis_cases;
    uint64_t high_recovered_symbols;
    uint64_t missing_cn_detected;
    uint64_t profiles;
    uint64_t mds_subsets;

    Counts()
        : field_checks(0)
        , subspace_checks(0)
        , normalization_checks(0)
        , high_basis_cases(0)
        , high_recovered_symbols(0)
        , missing_cn_detected(0)
        , profiles(0)
        , mds_subsets(0)
    {}
};

void require(bool condition, const std::string& message)
{
    if (!condition)
        throw std::runtime_error(message);
}

std::string context(unsigned n, unsigned t, uint32_t mask, unsigned basis_index)
{
    std::ostringstream stream;
    stream << "n=" << n << " t=" << t << " erasures=0x" << std::hex << mask
           << std::dec << " basis=" << basis_index;
    return stream.str();
}

unsigned popcount(uint32_t value)
{
    unsigned result = 0;
    while (value != 0)
    {
        result += value & 1u;
        value >>= 1;
    }
    return result;
}

void test_field(const BinaryField& field, bool exhaustive_products, Counts* counts)
{
    const unsigned order = field.order();
    for (unsigned a = 0; a < order; ++a)
    {
        require(field.add(static_cast<Element>(a), static_cast<Element>(a)) == 0,
            "field characteristic is not two");
        require(field.multiply(static_cast<Element>(a), 0) == 0,
            "zero multiplication failed");
        require(field.multiply(static_cast<Element>(a), 1) == a,
            "multiplicative identity failed");
        require(
            field.polynomial_to_coordinate(field.coordinate_to_polynomial(static_cast<Element>(a))) == a,
            "coordinate basis round trip failed");
        counts->field_checks += 4;
        if (a != 0)
        {
            require(field.multiply(static_cast<Element>(a), field.inverse(static_cast<Element>(a))) == 1,
                "field inverse failed");
            ++counts->field_checks;
        }
    }

    const unsigned product_limit = exhaustive_products ? order : 256;
    for (unsigned a_index = 0; a_index < product_limit; ++a_index)
    {
        const Element a = exhaustive_products
            ? static_cast<Element>(a_index)
            : static_cast<Element>((a_index * 40503u + 17u) & (order - 1));
        for (unsigned b_index = 0; b_index < product_limit; ++b_index)
        {
            const Element b = exhaustive_products
                ? static_cast<Element>(b_index)
                : static_cast<Element>((b_index * 17389u + a_index * 257u + 29u) & (order - 1));
            const Element polynomial_product = field.multiply_polynomial_symbols(
                field.coordinate_to_polynomial(a),
                field.coordinate_to_polynomial(b));
            require(field.multiply(a, b) == field.polynomial_to_coordinate(polynomial_product),
                "coordinate multiplication disagrees with polynomial reduction");
            require(field.multiply(a, b) == field.multiply(b, a),
                "field multiplication is not commutative");
            require(
                field.coordinate_to_polynomial(field.add(a, b)) ==
                    static_cast<Element>(field.coordinate_to_polynomial(a) ^
                                         field.coordinate_to_polynomial(b)),
                "coordinate XOR does not represent field addition");
            counts->field_checks += 3;
        }
    }
}

void test_gf4_field_laws(const BinaryField& field, Counts* counts)
{
    for (unsigned a = 0; a < field.order(); ++a)
        for (unsigned b = 0; b < field.order(); ++b)
            for (unsigned c = 0; c < field.order(); ++c)
            {
                const Element left = field.multiply(
                    static_cast<Element>(a),
                    field.add(static_cast<Element>(b), static_cast<Element>(c)));
                const Element right = field.add(
                    field.multiply(static_cast<Element>(a), static_cast<Element>(b)),
                    field.multiply(static_cast<Element>(a), static_cast<Element>(c)));
                require(left == right, "GF(2^4) distributive law failed");
                ++counts->field_checks;
            }
}

void test_subspaces(const BinaryField& field, Counts* counts)
{
    for (unsigned n = 0; n <= field.bits(); ++n)
    {
        const unsigned active_size = 1u << n;
        const Polynomial vanishing = subspace_vanishing_polynomial(field, n);
        const Polynomial derivative = polynomial_derivative(vanishing);
        const Element expected_derivative = subspace_derivative_constant(field, n);
        require(vanishing.size() == active_size + 1 && vanishing.back() == 1,
            "subspace polynomial has the wrong degree or leading coefficient");
        require(derivative.size() == 1 && derivative[0] == expected_derivative,
            "subspace derivative is not the nonzero-product constant c_n");
        require(expected_derivative != 0, "subspace derivative constant is zero");
        counts->subspace_checks += 3;

        for (unsigned x = 0; x < field.order(); ++x)
        {
            const bool is_root = polynomial_evaluate(field, vanishing, static_cast<Element>(x)) == 0;
            require(is_root == (x < active_size), "subspace polynomial has an incorrect root set");
            ++counts->subspace_checks;
        }
        for (unsigned x = 0; x < field.order(); ++x)
            for (unsigned y = 0; y < field.order(); ++y)
            {
                const Element left = polynomial_evaluate(
                    field, vanishing, static_cast<Element>(x ^ y));
                const Element right = static_cast<Element>(
                    polynomial_evaluate(field, vanishing, static_cast<Element>(x)) ^
                    polynomial_evaluate(field, vanishing, static_cast<Element>(y)));
                require(left == right, "subspace polynomial is not additive");
                ++counts->subspace_checks;
            }
    }
}

void test_active_normalization_identity(const BinaryField& field, Counts* counts)
{
    // R10 (10), tested with every active V_n rather than only the full field:
    // s_n = p_(N-T) Xbar_(N-T) s_t
    //       + sum_{i=t}^{n-1} p_(N-2^i) Xbar_(N-2^i) s_i(v_i).
    for (unsigned n = 1; n <= field.bits(); ++n)
    {
        const unsigned active_size = 1u << n;
        const LchBasis basis = make_lch_basis(field, n);
        for (unsigned t = 0; t < n; ++t)
        {
            const unsigned transform_size = 1u << t;
            const unsigned dimension = active_size - transform_size;
            Polynomial right = polynomial_scale(
                field,
                polynomial_multiply(
                    field,
                    basis.polynomials[dimension],
                    basis.subspace_polynomials[t]),
                basis.normalizers[dimension]);
            for (unsigned i = t; i < n; ++i)
            {
                const unsigned index = active_size - (1u << i);
                const Element subspace_at_basis = polynomial_evaluate(
                    field,
                    basis.subspace_polynomials[i],
                    static_cast<Element>(1u << i));
                const Element scale = field.multiply(
                    basis.normalizers[index], subspace_at_basis);
                right = polynomial_add(
                    right,
                    polynomial_scale(field, basis.polynomials[index], scale));
            }
            require(right == basis.subspace_polynomials[n],
                "active-parent normalization identity failed");
            ++counts->normalization_checks;
        }
    }
}

std::vector<Polynomial> cardinal_polynomials(
    const BinaryField& field,
    unsigned active_size)
{
    std::vector<Element> points(active_size, 0);
    for (unsigned i = 0; i < active_size; ++i)
        points[i] = static_cast<Element>(i);
    std::vector<Polynomial> cardinal(active_size);
    std::vector<Element> values(active_size, 0);
    for (unsigned i = 0; i < active_size; ++i)
    {
        values[i] = 1;
        cardinal[i] = lagrange_interpolate(field, points, values);
        values[i] = 0;
        for (unsigned x = 0; x < active_size; ++x)
            require(polynomial_evaluate(field, cardinal[i], static_cast<Element>(x)) == (x == i),
                "cardinal polynomial construction failed");
    }
    return cardinal;
}

Polynomial interpolate_masked_monomial(
    const BinaryField& field,
    const std::vector<Polynomial>& cardinal,
    uint32_t erasure_mask,
    unsigned monomial_degree)
{
    Polynomial result(1, 0);
    for (unsigned point = 0; point < cardinal.size(); ++point)
    {
        if ((erasure_mask & (1u << point)) != 0)
            continue;
        const Element value = field.power(static_cast<Element>(point), monomial_degree);
        result = polynomial_add(result, polynomial_scale(field, cardinal[point], value));
    }
    return result;
}

Polynomial make_locator(
    const BinaryField& field,
    unsigned active_size,
    uint32_t erasure_mask)
{
    Polynomial locator(1, 1);
    for (unsigned point = 0; point < active_size; ++point)
    {
        if ((erasure_mask & (1u << point)) == 0)
            continue;
        Polynomial factor(2, 0);
        factor[0] = static_cast<Element>(point);
        factor[1] = 1;
        locator = polynomial_multiply(field, locator, factor);
    }
    return locator;
}

Element locator_derivative_at(
    const BinaryField& field,
    unsigned active_size,
    uint32_t erasure_mask,
    unsigned erased_point)
{
    Element result = 1;
    for (unsigned point = 0; point < active_size; ++point)
    {
        if (point != erased_point && (erasure_mask & (1u << point)) != 0)
            result = field.multiply(result, static_cast<Element>(point ^ erased_point));
    }
    return result;
}

void test_active_high_congruence(const BinaryField& field, Counts* counts)
{
    // Every erasure set of the required size and every monomial basis vector is
    // enough to prove these linear identities for every degree-below-D f.
    for (unsigned n = 1; n <= field.bits(); ++n)
    {
        const unsigned active_size = 1u << n;
        const Polynomial vanishing = subspace_vanishing_polynomial(field, n);
        const Element active_derivative = subspace_derivative_constant(field, n);
        const LchBasis basis = make_lch_basis(field, n);
        const std::vector<Polynomial> cardinal = cardinal_polynomials(field, active_size);

        for (unsigned t = 0; t < n; ++t)
        {
            const unsigned transform_size = 1u << t;
            const unsigned dimension = active_size - transform_size;
            const Element p_dimension = basis.normalizers[dimension];
            const Element transform_derivative = subspace_derivative_constant(field, t);
            const uint32_t mask_limit = 1u << active_size;

            for (uint32_t erasure_mask = 0; erasure_mask < mask_limit; ++erasure_mask)
            {
                if (popcount(erasure_mask) != transform_size)
                    continue;
                const Polynomial locator = make_locator(field, active_size, erasure_mask);

                for (unsigned monomial_degree = 0;
                     monomial_degree < dimension;
                     ++monomial_degree)
                {
                    Polynomial original(monomial_degree + 1, 0);
                    original[monomial_degree] = 1;
                    const Polynomial masked = interpolate_masked_monomial(
                        field, cardinal, erasure_mask, monomial_degree);

                    // f Lambda = ftilde Lambda + q s_n on the active parent.
                    const Polynomial congruence_difference = polynomial_add(
                        polynomial_multiply(field, original, locator),
                        polynomial_multiply(field, masked, locator));
                    const Polynomial q = polynomial_divide_exact(
                        field, congruence_difference, vanishing);
                    require(q.size() <= transform_size,
                        "active high quotient has degree >= T: " +
                        context(n, t, erasure_mask, monomial_degree));

                    // Divide ftilde at the Xbar_D coefficient boundary.
                    const std::vector<Element> masked_lch = polynomial_to_lch_coefficients(
                        field, basis, masked);
                    std::vector<Element> h_coefficients(transform_size, 0);
                    for (unsigned i = 0; i < transform_size; ++i)
                        h_coefficients[i] = masked_lch[dimension + i];
                    const Polynomial h = lch_coefficients_to_polynomial(
                        field, basis, h_coefficients);

                    // The extracted key polynomial must have degree below T.
                    const Polynomial z = polynomial_add(
                        polynomial_multiply(field, h, locator),
                        polynomial_scale(
                            field,
                            polynomial_multiply(
                                field, q, basis.subspace_polynomials[t]),
                            p_dimension));
                    require(z.size() <= transform_size,
                        "active high key polynomial has degree >= T: " +
                        context(n, t, erasure_mask, monomial_degree));

                    for (unsigned point = 0; point < active_size; ++point)
                    {
                        if ((erasure_mask & (1u << point)) == 0)
                            continue;
                        const Element locator_derivative = locator_derivative_at(
                            field, active_size, erasure_mask, point);
                        require(locator_derivative != 0, "locator derivative is zero");
                        const Element expected = polynomial_evaluate(
                            field, original, static_cast<Element>(point));
                        Element recovered = 0;
                        Element recovered_without_cn = 0;

                        if (point >= transform_size)
                        {
                            const Element s_t = polynomial_evaluate(
                                field,
                                basis.subspace_polynomials[t],
                                static_cast<Element>(point));
                            const Element denominator = field.multiply(
                                field.multiply(p_dimension, s_t), locator_derivative);
                            recovered_without_cn = field.divide(
                                polynomial_evaluate(field, z, static_cast<Element>(point)),
                                denominator);
                            recovered = field.multiply(active_derivative, recovered_without_cn);
                        }
                        else
                        {
                            const Element z_derivative = polynomial_evaluate(
                                field,
                                polynomial_derivative(z),
                                static_cast<Element>(point));
                            const Element h_locator_derivative = field.multiply(
                                polynomial_evaluate(field, h, static_cast<Element>(point)),
                                locator_derivative);
                            const Element numerator = static_cast<Element>(
                                z_derivative ^ h_locator_derivative);
                            const Element denominator = field.multiply(
                                field.multiply(p_dimension, transform_derivative),
                                locator_derivative);
                            recovered_without_cn = field.divide(numerator, denominator);
                            recovered = field.multiply(active_derivative, recovered_without_cn);
                        }
                        require(recovered == expected,
                            "active high c_n recovery factor failed: " +
                            context(n, t, erasure_mask, monomial_degree));
                        if (active_derivative != 1 && recovered_without_cn != expected)
                            ++counts->missing_cn_detected;
                        ++counts->high_recovered_symbols;
                    }
                    ++counts->high_basis_cases;
                }
            }
        }
    }
    require(counts->missing_cn_detected != 0,
        "small-field suite did not detect omission of the active c_n factor");
}

void verify_inverse(
    const BinaryField& field,
    const Matrix& selected,
    const Matrix& inverse)
{
    const unsigned n = static_cast<unsigned>(selected.size());
    for (unsigned row = 0; row < n; ++row)
        for (unsigned column = 0; column < n; ++column)
        {
            Element product = 0;
            for (unsigned inner = 0; inner < n; ++inner)
                product ^= field.multiply(selected[row][inner], inverse[inner][column]);
            require(product == (row == column), "matrix inverse verification failed");
        }
}

void test_profile_mds(const BinaryField& field, ProfileKind kind, Counts* counts)
{
    for (unsigned k = 1; k < field.order(); ++k)
        for (unsigned r = 1; r < field.order(); ++r)
        {
            const ProfileLayout layout = make_profile_layout(kind, k, r);
            if (layout.parent_size > field.order())
                continue;
            const Matrix generator = direct_systematic_generator(field, layout);
            require(generator.size() == k + r, "generator has the wrong row count");
            for (unsigned row = 0; row < k; ++row)
                for (unsigned column = 0; column < k; ++column)
                    require(generator[row][column] == (row == column),
                        "generator is not systematic");

            // Cross-check matrix construction against direct interpolation of
            // every systematic basis word, including shortened zero positions.
            std::vector<Element> systematic_points(layout.parent_dimension, 0);
            for (unsigned i = 0; i < layout.parent_dimension; ++i)
                systematic_points[i] = static_cast<Element>(layout.systematic_coordinates[i]);
            std::vector<Element> values(layout.parent_dimension, 0);
            for (unsigned column = 0; column < k; ++column)
            {
                values[column] = 1;
                const Polynomial interpolated = lagrange_interpolate(field, systematic_points, values);
                for (unsigned parity = 0; parity < r; ++parity)
                    require(
                        generator[k + parity][column] == polynomial_evaluate(
                            field,
                            interpolated,
                            static_cast<Element>(layout.parity_coordinates[parity])),
                        "generator disagrees with polynomial interpolation");
                values[column] = 0;
            }

            std::vector<Element> message(k, 0);
            for (unsigned i = 0; i < k; ++i)
                message[i] = static_cast<Element>(i + 1);
            const std::vector<Element> encoded = matrix_vector_multiply(field, generator, message);
            const uint32_t subset_limit = 1u << (k + r);
            for (uint32_t subset = 0; subset < subset_limit; ++subset)
            {
                if (popcount(subset) != k)
                    continue;
                Matrix selected;
                std::vector<unsigned> rows;
                std::vector<Element> received;
                for (unsigned row = 0; row < k + r; ++row)
                {
                    if ((subset & (1u << row)) == 0)
                        continue;
                    selected.push_back(generator[row]);
                    rows.push_back(row);
                    received.push_back(encoded[row]);
                }
                Matrix inverse;
                require(invert_matrix(field, selected, &inverse),
                    "transmitted coordinate subset is not MDS");
                verify_inverse(field, selected, inverse);
                require(direct_recover(field, generator, rows, received) == message,
                    "direct recovery did not reproduce the message");
                ++counts->mds_subsets;
            }
            ++counts->profiles;
        }
}

} // namespace

int main()
{
    try
    {
        Counts counts;
        const BinaryField gf4 = make_gf4();
        const BinaryField gf8 = make_legacy_gf8();
        const BinaryField gf16 = make_legacy_gf16();

        test_field(gf4, true, &counts);
        test_field(gf8, true, &counts);
        test_field(gf16, false, &counts);
        test_gf4_field_laws(gf4, &counts);
        test_subspaces(gf4, &counts);
        test_active_normalization_identity(gf4, &counts);
        test_active_high_congruence(gf4, &counts);
        test_profile_mds(gf4, kLegacyHigh, &counts);
        test_profile_mds(gf4, kLow, &counts);

        std::cout
            << "PASS direct_oracle"
            << " field_checks=" << counts.field_checks
            << " subspace_checks=" << counts.subspace_checks
            << " normalization_checks=" << counts.normalization_checks
            << " high_basis_cases=" << counts.high_basis_cases
            << " high_recovered_symbols=" << counts.high_recovered_symbols
            << " missing_cn_detected=" << counts.missing_cn_detected
            << " profiles=" << counts.profiles
            << " mds_subsets=" << counts.mds_subsets
            << std::endl;
        return 0;
    }
    catch (const std::exception& exception)
    {
        std::cerr << "FAIL direct_oracle: " << exception.what() << std::endl;
        return 1;
    }
}
