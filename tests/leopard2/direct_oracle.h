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

#include <vector>

namespace leopard2_test {

typedef uint16_t Element;
typedef std::vector<Element> Polynomial; // Monomial basis, low coefficient first.
typedef std::vector<std::vector<Element> > Matrix;

// Slow, test-only binary-extension field.  Public Elements are coordinate bits
// in the supplied vector-space basis; multiplication converts to a polynomial
// basis and performs shift-and-reduce arithmetic without Leopard lookup tables.
class BinaryField
{
public:
    BinaryField(
        unsigned bits,
        uint32_t irreducible_polynomial,
        const std::vector<Element>& coordinate_basis);

    unsigned bits() const { return bits_; }
    unsigned order() const { return order_; }
    uint32_t polynomial() const { return polynomial_; }
    const std::vector<Element>& coordinate_basis() const { return basis_; }

    Element add(Element a, Element b) const { return static_cast<Element>(a ^ b); }
    Element multiply(Element a, Element b) const;
    Element power(Element a, uint32_t exponent) const;
    Element inverse(Element a) const;
    Element divide(Element a, Element b) const;

    Element coordinate_to_polynomial(Element value) const;
    Element polynomial_to_coordinate(Element value) const;
    Element multiply_polynomial_symbols(Element a, Element b) const;

private:
    unsigned bits_;
    unsigned order_;
    uint32_t polynomial_;
    std::vector<Element> basis_;
    std::vector<Element> coordinate_to_polynomial_;
    std::vector<Element> polynomial_to_coordinate_;
};

BinaryField make_gf4();
BinaryField make_legacy_gf8();
BinaryField make_legacy_gf16();

Polynomial polynomial_add(const Polynomial& a, const Polynomial& b);
Polynomial polynomial_scale(
    const BinaryField& field,
    const Polynomial& polynomial,
    Element scalar);
Polynomial polynomial_multiply(
    const BinaryField& field,
    const Polynomial& a,
    const Polynomial& b);
Polynomial polynomial_derivative(const Polynomial& polynomial);
Element polynomial_evaluate(
    const BinaryField& field,
    const Polynomial& polynomial,
    Element x);
Polynomial polynomial_divide_exact(
    const BinaryField& field,
    const Polynomial& numerator,
    const Polynomial& denominator);

Polynomial lagrange_interpolate(
    const BinaryField& field,
    const std::vector<Element>& points,
    const std::vector<Element>& values);
Element lagrange_evaluate(
    const BinaryField& field,
    const std::vector<Element>& points,
    const std::vector<Element>& values,
    Element x);

bool invert_matrix(
    const BinaryField& field,
    const Matrix& input,
    Matrix* inverse);
std::vector<Element> matrix_vector_multiply(
    const BinaryField& field,
    const Matrix& matrix,
    const std::vector<Element>& vector);

enum ProfileKind
{
    kLegacyHigh,
    kLow
};

struct ProfileLayout
{
    ProfileKind kind;
    unsigned original_count;
    unsigned recovery_count;
    unsigned parent_size;
    unsigned parent_dimension;
    unsigned padded_side; // T for high, P for low.
    std::vector<unsigned> systematic_coordinates; // Includes shortened zeros.
    std::vector<unsigned> parity_coordinates; // Transmitted parity only.
};

ProfileLayout make_profile_layout(
    ProfileKind kind,
    unsigned original_count,
    unsigned recovery_count);

// Public row order is K systematic symbols followed by R transmitted parities.
Matrix direct_systematic_generator(
    const BinaryField& field,
    const ProfileLayout& layout);
std::vector<Element> direct_encode(
    const BinaryField& field,
    const ProfileLayout& layout,
    const std::vector<Element>& message);
std::vector<Element> direct_recover(
    const BinaryField& field,
    const Matrix& generator,
    const std::vector<unsigned>& received_rows,
    const std::vector<Element>& received_values);

// Independent polynomial construction of V_n's vanishing polynomial.
Polynomial subspace_vanishing_polynomial(
    const BinaryField& field,
    unsigned dimension);
Element subspace_derivative_constant(
    const BinaryField& field,
    unsigned dimension);

// Normalized Lin-Chung-Han basis polynomials Xbar_i and p_i values for V_n.
struct LchBasis
{
    std::vector<Polynomial> subspace_polynomials;
    std::vector<Element> normalizers;
    std::vector<Polynomial> polynomials;
};

LchBasis make_lch_basis(const BinaryField& field, unsigned active_dimension);
std::vector<Element> polynomial_to_lch_coefficients(
    const BinaryField& field,
    const LchBasis& basis,
    const Polynomial& polynomial);
Polynomial lch_coefficients_to_polynomial(
    const BinaryField& field,
    const LchBasis& basis,
    const std::vector<Element>& coefficients);

} // namespace leopard2_test
