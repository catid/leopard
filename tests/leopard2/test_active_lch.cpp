/*
    Copyright (c) 2017 Christopher A. Taylor.  All rights reserved.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the conditions in the Leopard-RS
    BSD license are met.  See the repository LICENSE file.
*/

#include "direct_oracle.h"
#include "LeopardFF8.h"
#include "LeopardFF16.h"
#include "leopard2.h"

#include <algorithm>
#include <stdint.h>
#include <cstring>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <utility>
#include <vector>

namespace {

using leopard2_test::BinaryField;
using leopard2_test::Element;
using leopard2_test::LchBasis;
using leopard2_test::Polynomial;

static const uint64_t kKernelBytes = 64;

struct Counts
{
    uint64_t skew_factors;
    uint64_t normalizers;
    uint64_t subspace_factors;
    uint64_t prepared_factors;
    uint64_t forward_symbols;
    uint64_t inverse_symbols;
    uint64_t independent_inverse_symbols;
    uint64_t derivative_symbols;
    uint64_t lane_symbols;
    uint64_t transforms;

    Counts()
        : skew_factors(0), normalizers(0), subspace_factors(0)
        , prepared_factors(0), forward_symbols(0), inverse_symbols(0)
        , independent_inverse_symbols(0), derivative_symbols(0)
        , lane_symbols(0), transforms(0)
    {}
};

void require(bool condition, const std::string& message)
{
    if (!condition)
        throw std::runtime_error(message);
}

void require_result(leo2_result result, const char* operation)
{
    if (result != LEO2_SUCCESS)
    {
        std::ostringstream stream;
        stream << operation << ": " << leo2_result_string(result);
        throw std::runtime_error(stream.str());
    }
}

unsigned log2_exact(unsigned value)
{
    unsigned result = 0;
    while ((1u << result) < value)
        ++result;
    require((1u << result) == value, "non-dyadic test size");
    return result;
}

Element direct_subspace_at(
    const BinaryField& field,
    unsigned size,
    Element point)
{
    Element result = 1;
    for (unsigned i = 0; i < size; ++i)
        result = field.multiply(result, static_cast<Element>(point ^ i));
    return result;
}

Element direct_subspace_derivative(
    const BinaryField& field,
    unsigned size)
{
    Element result = 1;
    for (unsigned i = 1; i < size; ++i)
        result = field.multiply(result, static_cast<Element>(i));
    return result;
}

std::vector<Element> direct_normalizer_factors(const BinaryField& field)
{
    std::vector<Element> factors(field.bits(), 1);
    for (unsigned bit = 0; bit < field.bits(); ++bit)
    {
        const unsigned width = 1u << bit;
        factors[bit] = direct_subspace_at(
            field, width, static_cast<Element>(width));
        require(factors[bit] != 0, "zero LCH basis normalizer");
    }
    return factors;
}

Element direct_lch_normalizer(
    const BinaryField& field,
    const std::vector<Element>& factors,
    unsigned index)
{
    Element result = 1;
    for (unsigned bit = 0; bit < field.bits(); ++bit)
        if ((index & (1u << bit)) != 0)
            result = field.multiply(result, factors[bit]);
    return result;
}

Element direct_lch_basis_value(
    const BinaryField& field,
    const std::vector<Element>& factors,
    unsigned index,
    Element point)
{
    Element result = 1;
    for (unsigned bit = 0; bit < field.bits(); ++bit)
    {
        if ((index & (1u << bit)) == 0)
            continue;
        const unsigned width = 1u << bit;
        const Element numerator = direct_subspace_at(field, width, point);
        result = field.multiply(result, field.divide(numerator, factors[bit]));
    }
    return result;
}

Element direct_sparse_lch_evaluate(
    const BinaryField& field,
    const std::vector<Element>& factors,
    const std::vector<Element>& coefficients,
    Element point)
{
    Element result = 0;
    for (unsigned i = 0; i < coefficients.size(); ++i)
    {
        if (coefficients[i] == 0)
            continue;
        result ^= field.multiply(coefficients[i],
            direct_lch_basis_value(field, factors, i, point));
    }
    return result;
}

typedef std::pair<unsigned, Element> SparseTerm;

std::vector<SparseTerm> sparse_terms(
    const std::vector<Element>& coefficients)
{
    std::vector<SparseTerm> result;
    for (unsigned i = 0; i < coefficients.size(); ++i)
        if (coefficients[i] != 0)
            result.push_back(SparseTerm(i, coefficients[i]));
    return result;
}

std::vector<Element> normalizer_inverses(
    const BinaryField& field,
    const std::vector<Element>& normalizer_factors)
{
    std::vector<Element> result(normalizer_factors.size(), 0);
    for (unsigned i = 0; i < normalizer_factors.size(); ++i)
        result[i] = field.inverse(normalizer_factors[i]);
    return result;
}

Element direct_sparse_lch_evaluate_recurrence(
    const BinaryField& field,
    const std::vector<Element>& normalizer_factors,
    const std::vector<Element>& inverse_normalizers,
    const std::vector<SparseTerm>& terms,
    Element point)
{
    Element normalized_subspaces[16];
    require(field.bits() <= 16, "test field exceeds recurrence storage");

    // V_(j+1) is V_j union (v_j + V_j), and s_j is additive.  Hence
    // s_(j+1)(x) = s_j(x) * (s_j(x) + s_j(v_j)).  The independently
    // product-generated inverse normalizers turn these into Xbar factors.
    Element subspace = point;
    for (unsigned bit = 0; bit < field.bits(); ++bit)
    {
        normalized_subspaces[bit] = field.multiply(
            subspace, inverse_normalizers[bit]);
        subspace = field.multiply(
            subspace, subspace ^ normalizer_factors[bit]);
    }

    Element result = 0;
    for (unsigned term_i = 0; term_i < terms.size(); ++term_i)
    {
        Element basis = 1;
        const unsigned index = terms[term_i].first;
        for (unsigned bit = 0; bit < field.bits(); ++bit)
            if ((index & (1u << bit)) != 0)
                basis = field.multiply(basis, normalized_subspaces[bit]);
        result ^= field.multiply(terms[term_i].second, basis);
    }
    return result;
}

std::vector<Element> direct_derivative_accumulation(
    const BinaryField& field,
    const std::vector<Element>& normalizer_factors,
    const std::vector<Element>& coefficients)
{
    std::vector<Element> result = coefficients;
    for (unsigned input = 1; input < coefficients.size(); ++input)
    {
        const Element coefficient = coefficients[input];
        if (coefficient == 0)
            continue;
        for (unsigned bit = 0; bit < field.bits(); ++bit)
        {
            const unsigned width = 1u << bit;
            if (width >= coefficients.size() || (input & width) == 0)
                continue;
            const Element derivative = direct_subspace_derivative(field, width);
            const Element scale = field.divide(
                derivative, normalizer_factors[bit]);
            result[input ^ width] ^= field.multiply(coefficient, scale);
        }
    }
    return result;
}

std::vector<Element> deterministic_coefficients(
    const BinaryField& field,
    unsigned size,
    uint64_t salt)
{
    std::vector<Element> result(size, 0);
    const unsigned mask = field.order() - 1;
    const unsigned positions[] = {
        0, 1, size / 4, size / 2, size > 1 ? size - 2 : 0, size - 1
    };
    for (unsigned i = 0; i < sizeof(positions) / sizeof(positions[0]); ++i)
    {
        const unsigned position = positions[i];
        salt += UINT64_C(0x9e3779b97f4a7c15);
        uint64_t value = salt;
        value = (value ^ (value >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
        value = (value ^ (value >> 27)) * UINT64_C(0x94d049bb133111eb);
        value ^= value >> 31;
        result[position] ^= static_cast<Element>((value & mask) | 1u);
    }
    return result;
}

class KernelWork
{
public:
    explicit KernelWork(unsigned size)
        : shards_(size, std::vector<uint8_t>(kKernelBytes, 0))
        , pointers_(size, NULL)
    {
        for (unsigned i = 0; i < size; ++i)
            pointers_[i] = &shards_[i][0];
    }

    std::vector<uint8_t>& shard(unsigned index) { return shards_[index]; }
    const std::vector<uint8_t>& shard(unsigned index) const { return shards_[index]; }
    void** pointers() { return &pointers_[0]; }

private:
    std::vector<std::vector<uint8_t> > shards_;
    std::vector<void*> pointers_;
};

struct Production8
{
    typedef leopard::ff8::ffe_t Log;
    static unsigned order() { return leopard::ff8::kOrder; }
    static void set(std::vector<uint8_t>& shard, Element value)
    {
        std::fill(shard.begin(), shard.end(), static_cast<uint8_t>(value));
    }
    static unsigned lane_count() { return 64; }
    static void set_lane(std::vector<uint8_t>& shard, unsigned lane, Element value)
    {
        shard[lane] = static_cast<uint8_t>(value);
    }
    static Element get_lane(const std::vector<uint8_t>& shard, unsigned lane)
    {
        return shard[lane];
    }
    static Element get(const std::vector<uint8_t>& shard)
    {
        const Element value = shard[0];
        for (unsigned i = 1; i < shard.size(); ++i)
            require(shard[i] == value, "GF8 SIMD lanes differ");
        return value;
    }
    static Element skew(unsigned index) { return leopard::ff8::TestOnlyFFTMultiplier(index); }
    static Element derivative(unsigned size) { return leopard::ff8::TestOnlySubspaceDerivative(size); }
    static Element subspace(unsigned size, unsigned shift) { return leopard::ff8::TestOnlySubspaceAt(size, shift); }
    static Element normalizer(unsigned index) { return leopard::ff8::TestOnlyLchNormalizer(index); }
    static Log log(Element value) { return leopard::ff8::ElementLog(static_cast<Log>(value)); }
    static void prepare_low(unsigned n, unsigned p, Log* output) { leopard::ff8::PrepareLowDecode(n, p, output); }
    static void prepare_high(unsigned n, unsigned t, Log* output) { leopard::ff8::PrepareHighDecode(n, t, output); }
    static void forward(unsigned n, unsigned shift, unsigned count, void** work)
    {
        leopard::ff8::TestOnlyLchForward(kKernelBytes, n, shift, count, work);
    }
    static void inverse(unsigned n, unsigned shift, unsigned count, void** work)
    {
        leopard::ff8::TestOnlyLchInverse(kKernelBytes, n, shift, count, work);
    }
    static void derivative_accumulate(unsigned n, void** work)
    {
        leopard::ff8::TestOnlyAddFormalDerivative(kKernelBytes, n, work);
    }
};

struct Production16
{
    typedef leopard::ff16::ffe_t Log;
    static unsigned order() { return leopard::ff16::kOrder; }
    static void set(std::vector<uint8_t>& shard, Element value)
    {
        for (unsigned lane = 0; lane < 32; ++lane)
        {
            shard[lane] = static_cast<uint8_t>(value);
            shard[32 + lane] = static_cast<uint8_t>(value >> 8);
        }
    }
    static unsigned lane_count() { return 32; }
    static void set_lane(std::vector<uint8_t>& shard, unsigned lane, Element value)
    {
        shard[lane] = static_cast<uint8_t>(value);
        shard[32 + lane] = static_cast<uint8_t>(value >> 8);
    }
    static Element get_lane(const std::vector<uint8_t>& shard, unsigned lane)
    {
        return static_cast<Element>(shard[lane] |
            (static_cast<unsigned>(shard[32 + lane]) << 8));
    }
    static Element get(const std::vector<uint8_t>& shard)
    {
        const Element value = static_cast<Element>(
            shard[0] | (static_cast<unsigned>(shard[32]) << 8));
        for (unsigned lane = 1; lane < 32; ++lane)
            require(static_cast<Element>(shard[lane] |
                    (static_cast<unsigned>(shard[32 + lane]) << 8)) == value,
                "GF16 SIMD lanes differ");
        return value;
    }
    static Element skew(unsigned index) { return leopard::ff16::TestOnlyFFTMultiplier(index); }
    static Element derivative(unsigned size) { return leopard::ff16::TestOnlySubspaceDerivative(size); }
    static Element subspace(unsigned size, unsigned shift) { return leopard::ff16::TestOnlySubspaceAt(size, shift); }
    static Element normalizer(unsigned index) { return leopard::ff16::TestOnlyLchNormalizer(index); }
    static Log log(Element value) { return leopard::ff16::ElementLog(static_cast<Log>(value)); }
    static void prepare_low(unsigned n, unsigned p, Log* output) { leopard::ff16::PrepareLowDecode(n, p, output); }
    static void prepare_high(unsigned n, unsigned t, Log* output) { leopard::ff16::PrepareHighDecode(n, t, output); }
    static void forward(unsigned n, unsigned shift, unsigned count, void** work)
    {
        leopard::ff16::TestOnlyLchForward(kKernelBytes, n, shift, count, work);
    }
    static void inverse(unsigned n, unsigned shift, unsigned count, void** work)
    {
        leopard::ff16::TestOnlyLchInverse(kKernelBytes, n, shift, count, work);
    }
    static void derivative_accumulate(unsigned n, void** work)
    {
        leopard::ff16::TestOnlyAddFormalDerivative(kKernelBytes, n, work);
    }
};

template<class Production>
void load_work(KernelWork* work, const std::vector<Element>& values)
{
    for (unsigned i = 0; i < values.size(); ++i)
        Production::set(work->shard(i), values[i]);
}

template<class Production>
std::vector<Element> read_work(const KernelWork& work, unsigned size)
{
    std::vector<Element> result(size, 0);
    for (unsigned i = 0; i < size; ++i)
        result[i] = Production::get(work.shard(i));
    return result;
}

template<class Production>
void test_constants(const BinaryField& field, Counts* counts)
{
    const unsigned order = Production::order();
    const std::vector<Element> normalizer_factors =
        direct_normalizer_factors(field);

    for (unsigned index = 0; index + 1 < order; ++index)
    {
        unsigned half = 1;
        while ((index & half) != 0)
            half <<= 1;
        const unsigned shift = index - (half - 1);
        const Element expected = field.divide(
            direct_subspace_at(field, half, static_cast<Element>(shift)),
            normalizer_factors[log2_exact(half)]);
        require(Production::skew(index) == expected,
            "production FFTSkew differs from direct subspace quotient");
        ++counts->skew_factors;
    }

    for (unsigned index = 0; index < order; ++index)
    {
        const Element actual = Production::normalizer(index);
        const Element expected = direct_lch_normalizer(
            field, normalizer_factors, index);
        if (actual != expected)
        {
            std::ostringstream stream;
            stream << "production LCH normalizer differs from direct product: GF"
                   << order << " index=" << index << " actual=" << actual
                   << " expected=" << expected;
            throw std::runtime_error(stream.str());
        }
        ++counts->normalizers;
    }

    for (unsigned n = 1; n <= order; n <<= 1)
    {
        require(Production::derivative(n) ==
                direct_subspace_derivative(field, n),
            "production subspace derivative differs from direct product");
        ++counts->subspace_factors;
        for (unsigned shift = 0; shift < order; shift += n)
        {
            require(Production::subspace(n, shift) ==
                    direct_subspace_at(field, n, static_cast<Element>(shift)),
                "production subspace evaluation differs from direct product");
            ++counts->subspace_factors;
        }

        for (unsigned side = 2; side < n; side <<= 1)
        {
            typedef typename Production::Log Log;
            std::vector<Log> low(n / side - 1, 0);
            Production::prepare_low(n, side, &low[0]);
            const Element numerator = direct_subspace_derivative(field, side);
            for (unsigned block = 1; block < n / side; ++block)
            {
                const unsigned shift = block * side;
                const Element expected = field.divide(numerator,
                    direct_subspace_at(field, side, static_cast<Element>(shift)));
                require(low[block - 1] == Production::log(expected),
                    "prepared low factor differs from direct c_k/s_k");
                ++counts->prepared_factors;
            }

            std::vector<Log> high(n, 0);
            Production::prepare_high(n, side, &high[0]);
            const Element denominator_normalizer = direct_lch_normalizer(
                field, normalizer_factors, n - side);
            const Element active_derivative = direct_subspace_derivative(field, n);
            for (unsigned coordinate = side; coordinate < n; coordinate += side)
            {
                const Element expected = field.divide(active_derivative,
                    field.multiply(denominator_normalizer,
                        direct_subspace_at(field, side,
                            static_cast<Element>(coordinate))));
                const Log expected_log = Production::log(expected);
                for (unsigned i = 0; i < side; ++i)
                {
                    require(high[coordinate + i] == expected_log,
                        "prepared high output factor differs from direct active-parent factor");
                    ++counts->prepared_factors;
                }
            }
        }

        if (n > order / 2)
            break;
    }
}

std::vector<Element> direct_inverse(
    const BinaryField& field,
    unsigned size,
    unsigned shift,
    const std::vector<Element>& values)
{
    std::vector<Element> points(size, 0);
    for (unsigned i = 0; i < size; ++i)
        points[i] = static_cast<Element>(shift ^ i);
    const Polynomial polynomial = leopard2_test::lagrange_interpolate(
        field, points, values);
    const LchBasis basis = leopard2_test::make_lch_basis(
        field, log2_exact(size));
    return leopard2_test::polynomial_to_lch_coefficients(
        field, basis, polynomial);
}

template<class Production>
void test_transform_case(
    const BinaryField& field,
    const std::vector<Element>& normalizer_factors,
    unsigned size,
    unsigned shift,
    unsigned direct_output_limit,
    bool direct_inverse_enabled,
    bool recurrence_inverse_enabled,
    Counts* counts)
{
    const std::vector<Element> coefficients = deterministic_coefficients(
        field, size, UINT64_C(0x4c43484143544956) ^
        (static_cast<uint64_t>(size) << 20) ^ shift);

    const unsigned requested_counts[] = {
        1, std::min<unsigned>(3, size),
        std::min<unsigned>(direct_output_limit, size)
    };
    for (unsigned count_i = 0;
         count_i < sizeof(requested_counts) / sizeof(requested_counts[0]);
         ++count_i)
    {
        const unsigned requested = requested_counts[count_i];
        if (count_i != 0 && requested == requested_counts[count_i - 1])
            continue;
        KernelWork work(size);
        load_work<Production>(&work, coefficients);
        Production::forward(size, shift, requested, work.pointers());
        for (unsigned output = 0; output < requested; ++output)
        {
            const Element expected = direct_sparse_lch_evaluate(
                field, normalizer_factors, coefficients,
                static_cast<Element>(shift ^ output));
            require(Production::get(work.shard(output)) == expected,
                "bare production LCH forward differs from direct evaluation");
            ++counts->forward_symbols;
        }
        ++counts->transforms;
    }

    KernelWork round_trip(size);
    load_work<Production>(&round_trip, coefficients);
    Production::forward(size, shift, size, round_trip.pointers());
    const std::vector<Element> values = read_work<Production>(round_trip, size);

    if (recurrence_inverse_enabled && !direct_inverse_enabled)
    {
        const unsigned strategic_outputs[] = { 0, size / 2, size - 1 };
        for (unsigned i = 0;
             i < sizeof(strategic_outputs) / sizeof(strategic_outputs[0]); ++i)
        {
            const unsigned output = strategic_outputs[i];
            const Element expected = direct_sparse_lch_evaluate(
                field, normalizer_factors, coefficients,
                static_cast<Element>(shift ^ output));
            require(values[output] == expected,
                "large production forward differs at strategic output");
            ++counts->forward_symbols;
        }
    }

    Production::inverse(size, shift, size, round_trip.pointers());
    const std::vector<Element> recovered = read_work<Production>(round_trip, size);
    require(recovered == coefficients, "bare LCH forward/inverse round trip failed");
    counts->inverse_symbols += size;
    counts->transforms += 2;

    if (direct_inverse_enabled)
    {
        const unsigned known_counts[] = {
            0, 1, size / 2, size > 1 ? size - 1 : 0, size
        };
        for (unsigned count_i = 0;
             count_i < sizeof(known_counts) / sizeof(known_counts[0]);
             ++count_i)
        {
            const unsigned known = known_counts[count_i];
            if (count_i != 0 && known == known_counts[count_i - 1])
                continue;
            std::vector<Element> prefix_values(size, 0);
            for (unsigned i = 0; i < known; ++i)
                prefix_values[i] = static_cast<Element>(
                    (i * 29u + size * 7u + shift * 3u + 1u) &
                    (field.order() - 1));
            const std::vector<Element> expected = direct_inverse(
                field, size, shift, prefix_values);
            KernelWork inverse_work(size);
            load_work<Production>(&inverse_work, prefix_values);
            Production::inverse(size, shift, known, inverse_work.pointers());
            require(read_work<Production>(inverse_work, size) == expected,
                "bare production LCH inverse differs from direct interpolation");
            counts->inverse_symbols += size;
            ++counts->transforms;
        }
    }
    else if (recurrence_inverse_enabled)
    {
        const std::vector<Element> inverse_normalizers =
            normalizer_inverses(field, normalizer_factors);
        const std::vector<SparseTerm> terms = sparse_terms(coefficients);
        std::vector<Element> direct_values(size, 0);
        for (unsigned output = 0; output < size; ++output)
            direct_values[output] = direct_sparse_lch_evaluate_recurrence(
                field, normalizer_factors, inverse_normalizers, terms,
                static_cast<Element>(shift ^ output));

        const unsigned strategic_outputs[] = { 0, size / 2, size - 1 };
        for (unsigned i = 0;
             i < sizeof(strategic_outputs) / sizeof(strategic_outputs[0]); ++i)
        {
            const unsigned output = strategic_outputs[i];
            require(direct_values[output] == direct_sparse_lch_evaluate(
                    field, normalizer_factors, coefficients,
                    static_cast<Element>(shift ^ output)),
                "subspace recurrence differs from direct product");
        }

        KernelWork independent_inverse(size);
        load_work<Production>(&independent_inverse, direct_values);
        Production::inverse(
            size, shift, size, independent_inverse.pointers());
        require(read_work<Production>(independent_inverse, size) == coefficients,
            "large LCH inverse differs from independent recurrence values");
        counts->inverse_symbols += size;
        counts->independent_inverse_symbols += size;
        ++counts->transforms;
    }

    const std::vector<Element> expected_derivative =
        direct_derivative_accumulation(field, normalizer_factors, coefficients);
    KernelWork derivative(size);
    load_work<Production>(&derivative, coefficients);
    Production::derivative_accumulate(size, derivative.pointers());
    require(read_work<Production>(derivative, size) == expected_derivative,
        "normalized-LCH f+f' accumulation differs from direct algebra");
    counts->derivative_symbols += size;
}

template<class Production>
void test_transforms(
    const BinaryField& field,
    bool every_shift,
    bool recurrence_inverse_enabled,
    Counts* counts)
{
    const std::vector<Element> normalizer_factors =
        direct_normalizer_factors(field);
    const unsigned order = Production::order();
    for (unsigned size = 2; size <= order; size <<= 1)
    {
        std::vector<unsigned> shifts;
        if (every_shift)
        {
            for (unsigned shift = 0; shift < order; shift += size)
                shifts.push_back(shift);
        }
        else
        {
            shifts.push_back(0);
            if (size < order)
                shifts.push_back(order - size);
        }
        for (unsigned i = 0; i < shifts.size(); ++i)
            test_transform_case<Production>(field, normalizer_factors,
                size, shifts[i], every_shift ? size : std::min<unsigned>(size, 3),
                size <= 16, recurrence_inverse_enabled, counts);
        if (size > order / 2)
            break;
    }
}

template<class Production>
void test_lane_varying_transform(
    const BinaryField& field,
    unsigned size,
    unsigned shift,
    Counts* counts)
{
    require(size >= 2 && size <= Production::order(),
        "invalid lane-varying transform size");
    require((shift & (size - 1)) == 0 && shift + size <= Production::order(),
        "invalid lane-varying transform shift");

    const std::vector<Element> factors = direct_normalizer_factors(field);
    const std::vector<Element> inverses = normalizer_inverses(field, factors);
    const unsigned lanes = Production::lane_count();
    std::vector<std::vector<Element> > coefficients(
        lanes, std::vector<Element>(size, 0));
    std::vector<std::vector<SparseTerm> > terms(lanes);
    KernelWork work(size);

    for (unsigned lane = 0; lane < lanes; ++lane)
    {
        coefficients[lane] = deterministic_coefficients(
            field, size, UINT64_C(0x4c414e4556415259) ^
            (static_cast<uint64_t>(size) << 24) ^
            (static_cast<uint64_t>(shift) << 7) ^ lane);
        terms[lane] = sparse_terms(coefficients[lane]);
        for (unsigned coefficient = 0; coefficient < size; ++coefficient)
            Production::set_lane(work.shard(coefficient), lane,
                coefficients[lane][coefficient]);
    }

    Production::forward(size, shift, size, work.pointers());
    for (unsigned lane = 0; lane < lanes; ++lane)
    {
        for (unsigned output = 0; output < size; ++output)
        {
            const Element expected = direct_sparse_lch_evaluate_recurrence(
                field, factors, inverses, terms[lane],
                static_cast<Element>(shift ^ output));
            require(Production::get_lane(work.shard(output), lane) == expected,
                "lane-varying LCH forward differs from direct recurrence");
            ++counts->forward_symbols;
            ++counts->lane_symbols;
        }

        const unsigned strategic_outputs[] = { 0, size / 2, size - 1 };
        for (unsigned i = 0;
             i < sizeof(strategic_outputs) / sizeof(strategic_outputs[0]); ++i)
        {
            const unsigned output = strategic_outputs[i];
            require(Production::get_lane(work.shard(output), lane) ==
                    direct_sparse_lch_evaluate(field, factors,
                        coefficients[lane],
                        static_cast<Element>(shift ^ output)),
                "lane recurrence differs from direct product");
        }
    }

    Production::inverse(size, shift, size, work.pointers());
    for (unsigned lane = 0; lane < lanes; ++lane)
        for (unsigned coefficient = 0; coefficient < size; ++coefficient)
        {
            require(Production::get_lane(work.shard(coefficient), lane) ==
                    coefficients[lane][coefficient],
                "lane-varying LCH inverse did not recover coefficients");
            ++counts->inverse_symbols;
            ++counts->lane_symbols;
        }
    counts->transforms += 2;
}

} // namespace

int main()
{
    try
    {
        leo2_context* context = NULL;
        require_result(leo2_context_create(NULL, &context), "context create");

        Counts counts;
        const BinaryField gf8 = leopard2_test::make_legacy_gf8();
        const BinaryField gf16 = leopard2_test::make_legacy_gf16();
        test_constants<Production8>(gf8, &counts);
        test_constants<Production16>(gf16, &counts);
        test_transforms<Production8>(gf8, true, false, &counts);
        test_transforms<Production16>(gf16, false, true, &counts);
        test_lane_varying_transform<Production8>(gf8, 64, 192, &counts);
        test_lane_varying_transform<Production16>(gf16, 4096, 0, &counts);
        test_lane_varying_transform<Production16>(
            gf16, 4096, Production16::order() - 4096, &counts);

        leo2_context_destroy(context);
        std::cout << "PASS active_lch"
                  << " skew_factors=" << counts.skew_factors
                  << " normalizers=" << counts.normalizers
                  << " subspace_factors=" << counts.subspace_factors
                  << " prepared_factors=" << counts.prepared_factors
                  << " forward_symbols=" << counts.forward_symbols
                  << " inverse_symbols=" << counts.inverse_symbols
                  << " independent_inverse_symbols="
                  << counts.independent_inverse_symbols
                  << " derivative_symbols=" << counts.derivative_symbols
                  << " lane_symbols=" << counts.lane_symbols
                  << " transforms=" << counts.transforms
                  << std::endl;
        return 0;
    }
    catch (const std::exception& error)
    {
        std::cerr << "FAIL active_lch: " << error.what() << std::endl;
        return 1;
    }
}
