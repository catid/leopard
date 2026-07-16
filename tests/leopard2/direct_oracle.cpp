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

#include <algorithm>
#include <limits>
#include <stdexcept>

namespace leopard2_test {
namespace {

void trim(Polynomial* polynomial)
{
    while (polynomial->size() > 1 && polynomial->back() == 0)
        polynomial->pop_back();
    if (polynomial->empty())
        polynomial->push_back(0);
}

bool is_zero(const Polynomial& polynomial)
{
    for (unsigned i = 0; i < polynomial.size(); ++i)
        if (polynomial[i] != 0)
            return false;
    return true;
}

unsigned next_power_of_two(unsigned value)
{
    unsigned result = 1;
    while (result < value)
    {
        if (result > std::numeric_limits<unsigned>::max() / 2)
            throw std::overflow_error("power-of-two overflow");
        result <<= 1;
    }
    return result;
}

Polynomial monomial_linear_factor(Element root)
{
    Polynomial result(2);
    result[0] = root; // x - root == x + root in characteristic two.
    result[1] = 1;
    return result;
}

} // namespace

BinaryField::BinaryField(
    unsigned bits,
    uint32_t irreducible_polynomial,
    const std::vector<Element>& coordinate_basis)
    : bits_(bits)
    , order_(bits < 31 ? (1u << bits) : 0)
    , polynomial_(irreducible_polynomial)
    , basis_(coordinate_basis)
{
    if (bits == 0 || bits > 16 || basis_.size() != bits || order_ == 0)
        throw std::invalid_argument("unsupported binary field size");
    if ((polynomial_ & (1u << bits_)) == 0)
        throw std::invalid_argument("field polynomial has the wrong degree");

    coordinate_to_polynomial_.assign(order_, 0);
    polynomial_to_coordinate_.assign(order_, std::numeric_limits<Element>::max());
    for (unsigned coordinate = 0; coordinate < order_; ++coordinate)
    {
        Element polynomial_value = 0;
        for (unsigned bit = 0; bit < bits_; ++bit)
            if ((coordinate & (1u << bit)) != 0)
                polynomial_value ^= basis_[bit];
        if (polynomial_value >= order_ ||
            polynomial_to_coordinate_[polynomial_value] != std::numeric_limits<Element>::max())
        {
            throw std::invalid_argument("coordinate basis is not linearly independent");
        }
        coordinate_to_polynomial_[coordinate] = polynomial_value;
        polynomial_to_coordinate_[polynomial_value] = static_cast<Element>(coordinate);
    }
}

Element BinaryField::coordinate_to_polynomial(Element value) const
{
    if (value >= order_)
        throw std::out_of_range("field coordinate out of range");
    return coordinate_to_polynomial_[value];
}

Element BinaryField::polynomial_to_coordinate(Element value) const
{
    if (value >= order_)
        throw std::out_of_range("polynomial-basis element out of range");
    return polynomial_to_coordinate_[value];
}

Element BinaryField::multiply_polynomial_symbols(Element a, Element b) const
{
    uint32_t product = 0;
    uint32_t left = a;
    uint32_t right = b;
    while (right != 0)
    {
        if ((right & 1u) != 0)
            product ^= left;
        right >>= 1;
        left <<= 1;
    }

    for (int bit = static_cast<int>(bits_ * 2 - 2);
         bit >= static_cast<int>(bits_);
         --bit)
    {
        if ((product & (1u << bit)) != 0)
            product ^= polynomial_ << (bit - bits_);
    }
    return static_cast<Element>(product);
}

Element BinaryField::multiply(Element a, Element b) const
{
    const Element polynomial_product = multiply_polynomial_symbols(
        coordinate_to_polynomial(a), coordinate_to_polynomial(b));
    return polynomial_to_coordinate(polynomial_product);
}

Element BinaryField::power(Element a, uint32_t exponent) const
{
    Element result = 1;
    while (exponent != 0)
    {
        if ((exponent & 1u) != 0)
            result = multiply(result, a);
        exponent >>= 1;
        if (exponent != 0)
            a = multiply(a, a);
    }
    return result;
}

Element BinaryField::inverse(Element a) const
{
    if (a == 0)
        throw std::domain_error("zero has no multiplicative inverse");
    return power(a, order_ - 2);
}

Element BinaryField::divide(Element a, Element b) const
{
    if (a == 0)
        return 0;
    return multiply(a, inverse(b));
}

BinaryField make_gf4()
{
    const Element basis[] = { 0x1, 0x2, 0x4, 0x8 };
    return BinaryField(4, 0x13, std::vector<Element>(basis, basis + 4));
}

BinaryField make_legacy_gf8()
{
    const Element basis[] = { 1, 214, 152, 146, 86, 200, 88, 230 };
    return BinaryField(8, 0x11d, std::vector<Element>(basis, basis + 8));
}

BinaryField make_legacy_gf16()
{
    const Element basis[] = {
        0x0001, 0xacca, 0x3c0e, 0x163e,
        0xc582, 0xed2e, 0x914c, 0x4012,
        0x6c98, 0x10d8, 0x6a72, 0xb900,
        0xfdb8, 0xfb34, 0xff38, 0x991e
    };
    return BinaryField(16, 0x1002d, std::vector<Element>(basis, basis + 16));
}

Polynomial polynomial_add(const Polynomial& a, const Polynomial& b)
{
    Polynomial result(std::max(a.size(), b.size()), 0);
    for (unsigned i = 0; i < a.size(); ++i)
        result[i] ^= a[i];
    for (unsigned i = 0; i < b.size(); ++i)
        result[i] ^= b[i];
    trim(&result);
    return result;
}

Polynomial polynomial_scale(
    const BinaryField& field,
    const Polynomial& polynomial,
    Element scalar)
{
    Polynomial result(polynomial.size(), 0);
    for (unsigned i = 0; i < polynomial.size(); ++i)
        result[i] = field.multiply(polynomial[i], scalar);
    trim(&result);
    return result;
}

Polynomial polynomial_multiply(
    const BinaryField& field,
    const Polynomial& a,
    const Polynomial& b)
{
    if (is_zero(a) || is_zero(b))
        return Polynomial(1, 0);
    Polynomial result(a.size() + b.size() - 1, 0);
    for (unsigned i = 0; i < a.size(); ++i)
        for (unsigned j = 0; j < b.size(); ++j)
            result[i + j] ^= field.multiply(a[i], b[j]);
    trim(&result);
    return result;
}

Polynomial polynomial_derivative(const Polynomial& polynomial)
{
    if (polynomial.size() <= 1)
        return Polynomial(1, 0);
    Polynomial result(polynomial.size() - 1, 0);
    // In characteristic two, only odd monomial powers survive.
    for (unsigned i = 1; i < polynomial.size(); i += 2)
        result[i - 1] = polynomial[i];
    trim(&result);
    return result;
}

Element polynomial_evaluate(
    const BinaryField& field,
    const Polynomial& polynomial,
    Element x)
{
    Element result = 0;
    for (unsigned i = static_cast<unsigned>(polynomial.size()); i-- > 0;)
        result = field.add(field.multiply(result, x), polynomial[i]);
    return result;
}

Polynomial polynomial_divide_exact(
    const BinaryField& field,
    const Polynomial& numerator,
    const Polynomial& denominator)
{
    Polynomial remainder = numerator;
    Polynomial divisor = denominator;
    trim(&remainder);
    trim(&divisor);
    if (is_zero(divisor))
        throw std::domain_error("polynomial division by zero");
    if (remainder.size() < divisor.size())
    {
        if (!is_zero(remainder))
            throw std::domain_error("polynomial division has a remainder");
        return Polynomial(1, 0);
    }

    Polynomial quotient(remainder.size() - divisor.size() + 1, 0);
    const Element inverse_leading = field.inverse(divisor.back());
    while (!is_zero(remainder) && remainder.size() >= divisor.size())
    {
        const unsigned shift = static_cast<unsigned>(remainder.size() - divisor.size());
        const Element factor = field.multiply(remainder.back(), inverse_leading);
        quotient[shift] ^= factor;
        for (unsigned i = 0; i < divisor.size(); ++i)
            remainder[i + shift] ^= field.multiply(factor, divisor[i]);
        trim(&remainder);
    }
    if (!is_zero(remainder))
        throw std::domain_error("polynomial division has a remainder");
    trim(&quotient);
    return quotient;
}

Polynomial lagrange_interpolate(
    const BinaryField& field,
    const std::vector<Element>& points,
    const std::vector<Element>& values)
{
    if (points.empty() || points.size() != values.size())
        throw std::invalid_argument("invalid interpolation vectors");
    Polynomial result(1, 0);
    for (unsigned i = 0; i < points.size(); ++i)
    {
        Polynomial cardinal(1, 1);
        Element denominator = 1;
        for (unsigned j = 0; j < points.size(); ++j)
        {
            if (i == j)
                continue;
            const Element difference = field.add(points[i], points[j]);
            if (difference == 0)
                throw std::invalid_argument("duplicate interpolation point");
            cardinal = polynomial_multiply(field, cardinal, monomial_linear_factor(points[j]));
            denominator = field.multiply(denominator, difference);
        }
        const Element scale = field.divide(values[i], denominator);
        result = polynomial_add(result, polynomial_scale(field, cardinal, scale));
    }
    trim(&result);
    return result;
}

Element lagrange_evaluate(
    const BinaryField& field,
    const std::vector<Element>& points,
    const std::vector<Element>& values,
    Element x)
{
    if (points.empty() || points.size() != values.size())
        throw std::invalid_argument("invalid interpolation vectors");
    Element result = 0;
    for (unsigned i = 0; i < points.size(); ++i)
    {
        Element numerator = values[i];
        Element denominator = 1;
        for (unsigned j = 0; j < points.size(); ++j)
        {
            if (i == j)
                continue;
            const Element difference = field.add(points[i], points[j]);
            if (difference == 0)
                throw std::invalid_argument("duplicate interpolation point");
            numerator = field.multiply(numerator, field.add(x, points[j]));
            denominator = field.multiply(denominator, difference);
        }
        result ^= field.divide(numerator, denominator);
    }
    return result;
}

bool invert_matrix(
    const BinaryField& field,
    const Matrix& input,
    Matrix* inverse)
{
    const unsigned n = static_cast<unsigned>(input.size());
    if (n == 0 || inverse == 0)
        throw std::invalid_argument("invalid matrix inversion request");
    Matrix augmented(n, std::vector<Element>(n * 2, 0));
    for (unsigned row = 0; row < n; ++row)
    {
        if (input[row].size() != n)
            throw std::invalid_argument("matrix is not square");
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
            return false;
        if (pivot != column)
            std::swap(augmented[pivot], augmented[column]);

        const Element scale = field.inverse(augmented[column][column]);
        for (unsigned j = 0; j < n * 2; ++j)
            augmented[column][j] = field.multiply(augmented[column][j], scale);

        for (unsigned row = 0; row < n; ++row)
        {
            if (row == column || augmented[row][column] == 0)
                continue;
            const Element factor = augmented[row][column];
            for (unsigned j = 0; j < n * 2; ++j)
                augmented[row][j] ^= field.multiply(factor, augmented[column][j]);
        }
    }

    inverse->assign(n, std::vector<Element>(n, 0));
    for (unsigned row = 0; row < n; ++row)
        for (unsigned column = 0; column < n; ++column)
            (*inverse)[row][column] = augmented[row][n + column];
    return true;
}

std::vector<Element> matrix_vector_multiply(
    const BinaryField& field,
    const Matrix& matrix,
    const std::vector<Element>& vector)
{
    std::vector<Element> result(matrix.size(), 0);
    for (unsigned row = 0; row < matrix.size(); ++row)
    {
        if (matrix[row].size() != vector.size())
            throw std::invalid_argument("matrix/vector dimensions do not match");
        for (unsigned column = 0; column < vector.size(); ++column)
            result[row] ^= field.multiply(matrix[row][column], vector[column]);
    }
    return result;
}

ProfileLayout make_profile_layout(
    ProfileKind kind,
    unsigned original_count,
    unsigned recovery_count)
{
    if (original_count == 0 || recovery_count == 0)
        throw std::invalid_argument("profile counts must be positive");

    ProfileLayout layout;
    layout.kind = kind;
    layout.original_count = original_count;
    layout.recovery_count = recovery_count;
    if (kind == kLegacyHigh)
    {
        const unsigned t = next_power_of_two(recovery_count);
        const unsigned n = next_power_of_two(original_count + t);
        layout.padded_side = t;
        layout.parent_size = n;
        layout.parent_dimension = n - t;
        for (unsigned i = t; i < n; ++i)
            layout.systematic_coordinates.push_back(i);
        for (unsigned i = 0; i < recovery_count; ++i)
            layout.parity_coordinates.push_back(i);
    }
    else
    {
        const unsigned p = next_power_of_two(original_count);
        const unsigned n = next_power_of_two(p + recovery_count);
        layout.padded_side = p;
        layout.parent_size = n;
        layout.parent_dimension = p;
        for (unsigned i = 0; i < p; ++i)
            layout.systematic_coordinates.push_back(i);
        for (unsigned i = 0; i < recovery_count; ++i)
            layout.parity_coordinates.push_back(p + i);
    }
    return layout;
}

Matrix direct_systematic_generator(
    const BinaryField& field,
    const ProfileLayout& layout)
{
    if (layout.parent_size > field.order() ||
        layout.systematic_coordinates.size() != layout.parent_dimension ||
        layout.parent_dimension < layout.original_count)
    {
        throw std::invalid_argument("profile does not fit the field");
    }

    Matrix generator(
        layout.original_count + layout.recovery_count,
        std::vector<Element>(layout.original_count, 0));
    for (unsigned i = 0; i < layout.original_count; ++i)
        generator[i][i] = 1;

    std::vector<Element> systematic_points(layout.parent_dimension, 0);
    for (unsigned i = 0; i < layout.parent_dimension; ++i)
        systematic_points[i] = static_cast<Element>(layout.systematic_coordinates[i]);

    std::vector<Element> values(layout.parent_dimension, 0);
    for (unsigned message_index = 0; message_index < layout.original_count; ++message_index)
    {
        values[message_index] = 1;
        for (unsigned parity_index = 0; parity_index < layout.recovery_count; ++parity_index)
        {
            generator[layout.original_count + parity_index][message_index] = lagrange_evaluate(
                field,
                systematic_points,
                values,
                static_cast<Element>(layout.parity_coordinates[parity_index]));
        }
        values[message_index] = 0;
    }
    return generator;
}

std::vector<Element> direct_encode(
    const BinaryField& field,
    const ProfileLayout& layout,
    const std::vector<Element>& message)
{
    if (message.size() != layout.original_count)
        throw std::invalid_argument("message length does not match profile");
    return matrix_vector_multiply(field, direct_systematic_generator(field, layout), message);
}

std::vector<Element> direct_recover(
    const BinaryField& field,
    const Matrix& generator,
    const std::vector<unsigned>& received_rows,
    const std::vector<Element>& received_values)
{
    if (generator.empty() || received_rows.size() != generator[0].size() ||
        received_rows.size() != received_values.size())
    {
        throw std::invalid_argument("recovery requires exactly K symbols");
    }
    Matrix selected(received_rows.size(), std::vector<Element>(received_rows.size(), 0));
    for (unsigned row = 0; row < received_rows.size(); ++row)
    {
        if (received_rows[row] >= generator.size())
            throw std::out_of_range("received generator row out of range");
        selected[row] = generator[received_rows[row]];
    }
    Matrix inverse;
    if (!invert_matrix(field, selected, &inverse))
        throw std::domain_error("selected code coordinates are singular");
    return matrix_vector_multiply(field, inverse, received_values);
}

Polynomial subspace_vanishing_polynomial(
    const BinaryField& field,
    unsigned dimension)
{
    if (dimension > field.bits())
        throw std::invalid_argument("subspace dimension exceeds field dimension");
    Polynomial result(1, 1);
    const unsigned count = 1u << dimension;
    for (unsigned i = 0; i < count; ++i)
        result = polynomial_multiply(field, result, monomial_linear_factor(static_cast<Element>(i)));
    return result;
}

Element subspace_derivative_constant(
    const BinaryField& field,
    unsigned dimension)
{
    if (dimension > field.bits())
        throw std::invalid_argument("subspace dimension exceeds field dimension");
    Element result = 1;
    const unsigned count = 1u << dimension;
    for (unsigned i = 1; i < count; ++i)
        result = field.multiply(result, static_cast<Element>(i));
    return result;
}

LchBasis make_lch_basis(const BinaryField& field, unsigned active_dimension)
{
    if (active_dimension > field.bits())
        throw std::invalid_argument("active LCH dimension exceeds field dimension");
    const unsigned n = 1u << active_dimension;
    LchBasis basis;
    basis.subspace_polynomials.resize(active_dimension + 1);
    for (unsigned i = 0; i <= active_dimension; ++i)
        basis.subspace_polynomials[i] = subspace_vanishing_polynomial(field, i);

    basis.normalizers.assign(n, 1);
    basis.polynomials.assign(n, Polynomial(1, 1));
    for (unsigned index = 0; index < n; ++index)
    {
        Element normalizer = 1;
        Polynomial polynomial(1, 1);
        for (unsigned bit = 0; bit < active_dimension; ++bit)
        {
            if ((index & (1u << bit)) == 0)
                continue;
            const Element factor = polynomial_evaluate(
                field,
                basis.subspace_polynomials[bit],
                static_cast<Element>(1u << bit));
            if (factor == 0)
                throw std::domain_error("invalid LCH coordinate basis");
            normalizer = field.multiply(normalizer, factor);
            polynomial = polynomial_multiply(field, polynomial, basis.subspace_polynomials[bit]);
        }
        basis.normalizers[index] = normalizer;
        basis.polynomials[index] = polynomial_scale(field, polynomial, field.inverse(normalizer));
    }
    return basis;
}

std::vector<Element> polynomial_to_lch_coefficients(
    const BinaryField& field,
    const LchBasis& basis,
    const Polynomial& polynomial)
{
    const unsigned n = static_cast<unsigned>(basis.polynomials.size());
    if (polynomial.size() > n)
        throw std::invalid_argument("polynomial exceeds LCH basis size");
    Polynomial remainder = polynomial;
    remainder.resize(n, 0);
    std::vector<Element> coefficients(n, 0);
    for (unsigned index = n; index-- > 0;)
    {
        const Polynomial& basis_polynomial = basis.polynomials[index];
        if (basis_polynomial.size() != index + 1 || basis_polynomial[index] == 0)
            throw std::domain_error("LCH basis is not triangular");
        const Element coefficient = field.divide(remainder[index], basis_polynomial[index]);
        coefficients[index] = coefficient;
        for (unsigned j = 0; j < basis_polynomial.size(); ++j)
            remainder[j] ^= field.multiply(coefficient, basis_polynomial[j]);
    }
    trim(&remainder);
    if (!is_zero(remainder))
        throw std::domain_error("LCH basis conversion left a remainder");
    return coefficients;
}

Polynomial lch_coefficients_to_polynomial(
    const BinaryField& field,
    const LchBasis& basis,
    const std::vector<Element>& coefficients)
{
    if (coefficients.size() > basis.polynomials.size())
        throw std::invalid_argument("too many LCH coefficients");
    Polynomial result(1, 0);
    for (unsigned i = 0; i < coefficients.size(); ++i)
        result = polynomial_add(
            result,
            polynomial_scale(field, basis.polynomials[i], coefficients[i]));
    return result;
}

} // namespace leopard2_test
