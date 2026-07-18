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

#include "LeopardFF8.h"
#include "LeopardFF16.h"

#include <algorithm>
#include <atomic>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <thread>
#include <vector>

namespace {

static void Require(bool condition, const char* message)
{
    if (!condition)
    {
        std::cerr << "locator test failed: " << message << std::endl;
        std::exit(1);
    }
}

static uint32_t Mix(uint32_t value)
{
    value ^= value >> 16;
    value *= 0x7feb352du;
    value ^= value >> 15;
    value *= 0x846ca68bu;
    return value ^ (value >> 16);
}

static unsigned Popcount32(uint32_t value)
{
    unsigned count = 0;
    while (value != 0)
    {
        value &= value - 1;
        ++count;
    }
    return count;
}

template<class Ffe>
static Ffe CanonicalLog(Ffe value)
{
    return value == std::numeric_limits<Ffe>::max() ? 0 : value;
}

template<class Ffe>
static bool EquivalentLogs(const std::vector<Ffe>& a, const std::vector<Ffe>& b)
{
    if (a.size() != b.size())
        return false;
    for (size_t i = 0; i < a.size(); ++i)
        if (CanonicalLog(a[i]) != CanonicalLog(b[i]))
            return false;
    return true;
}

static unsigned AddModulo(unsigned a, unsigned b, unsigned modulus);

template<class Ffe>
static uint64_t CheckSampledDirectCoordinates(
    unsigned n,
    const std::vector<uint8_t>& erasures,
    const std::vector<Ffe>& actual,
    uint32_t seed,
    Ffe (*element_log)(Ffe),
    unsigned modulus)
{
    const unsigned sample_limit = std::min(64u, n);
    std::vector<uint8_t> selected(n, 0);
    std::vector<unsigned> coordinates;
    coordinates.reserve(sample_limit);

    const unsigned required[] = { 0, n - 1, Mix(seed) & (n - 1) };
    for (unsigned i = 0;
         i < sizeof(required) / sizeof(required[0]) &&
             coordinates.size() < sample_limit;
         ++i)
    {
        if (!selected[required[i]])
        {
            selected[required[i]] = 1;
            coordinates.push_back(required[i]);
        }
    }
    const unsigned start = Mix(seed ^ 0xa511e9b3u) & (n - 1);
    const unsigned step = (Mix(seed ^ 0x63d83595u) | 1u) & (n - 1);
    for (unsigned attempt = 0; coordinates.size() < sample_limit; ++attempt)
    {
        const unsigned coordinate = (start + attempt * step) & (n - 1);
        if (!selected[coordinate])
        {
            selected[coordinate] = 1;
            coordinates.push_back(coordinate);
        }
    }

    bool crossed_modulus = false;
    for (size_t sample = 0; sample < coordinates.size(); ++sample)
    {
        const unsigned coordinate = coordinates[sample];
        unsigned expected = 0;
        uint64_t unreduced = 0;
        for (unsigned erased = 0; erased < n; ++erased)
        {
            if (!erasures[erased] || erased == coordinate)
                continue;
            const unsigned logarithm = static_cast<unsigned>(
                element_log(static_cast<Ffe>(coordinate ^ erased)));
            unreduced += logarithm;
            expected = AddModulo(expected, logarithm, modulus);
        }
        crossed_modulus = crossed_modulus || unreduced >= modulus;
        Require(CanonicalLog(actual[coordinate]) == expected,
            "sampled direct/active locator mismatch");
    }
    if (static_cast<uint64_t>(n) *
            std::count(erasures.begin(), erasures.end(), uint8_t(1)) >
        2u * 1024u * 1024u)
    {
        Require(crossed_modulus,
            "dense sampled locator did not exercise modular wrap");
    }
    return coordinates.size();
}

template<class Ffe, class Prepare, class PreparePermanent, class Direct,
    class Active, class Reference>
static uint64_t CheckPattern(
    unsigned n,
    unsigned erasure_count,
    uint32_t seed,
    Prepare prepare,
    PreparePermanent prepare_permanent,
    Direct direct,
    Active active,
    Reference reference,
    Ffe (*element_log)(Ffe),
    unsigned modulus)
{
    std::vector<uint8_t> erasures(n, 0);
    std::vector<unsigned> order(n);
    for (unsigned i = 0; i < n; ++i)
        order[i] = i;
    std::sort(order.begin(), order.end(), [seed](unsigned a, unsigned b) {
        const uint32_t ma = Mix(a ^ seed);
        const uint32_t mb = Mix(b ^ seed);
        return ma == mb ? a < b : ma < mb;
    });
    for (unsigned i = 0; i < erasure_count; ++i)
        erasures[order[i]] = 1;

    std::vector<Ffe> expected(n), direct_actual(n), active_actual(n), actual(n),
        base(n), with_base(n), with_dynamic_only(n);
    reference(n, &erasures[0], &expected[0]);
    const bool checked_direct =
        static_cast<uint64_t>(n) * erasure_count <= 2u * 1024u * 1024u;
    if (checked_direct)
    {
        direct(n, &erasures[0], &direct_actual[0]);
        Require(EquivalentLogs(direct_actual, expected),
            "direct/full locator mismatch");
    }
    active(n, &erasures[0], &active_actual[0]);
    Require(EquivalentLogs(active_actual, expected),
        "active/full locator mismatch");
    uint64_t sampled = 0;
    if (!checked_direct)
    {
        sampled = CheckSampledDirectCoordinates(
            n, erasures, active_actual, seed, element_log, modulus);
    }
    prepare(n, &erasures[0], &actual[0]);
    if (!EquivalentLogs(actual, expected))
    {
        size_t mismatch = 0;
        while (mismatch < actual.size() &&
               CanonicalLog(actual[mismatch]) == CanonicalLog(expected[mismatch]))
            ++mismatch;
        std::cerr << "direct mismatch n=" << n << " erasures=" << erasure_count
                  << " coordinate=" << mismatch << " expected="
                  << static_cast<unsigned>(expected[mismatch]) << " actual="
                  << static_cast<unsigned>(actual[mismatch]) << std::endl;
    }
    Require(EquivalentLogs(actual, expected), "sparse/full locator mismatch");

    std::vector<uint8_t> permanent(n, 0);
    unsigned selected = 0;
    for (unsigned i = 0; i < n; ++i)
    {
        if (erasures[i] && (selected++ % 3u) == 0)
            permanent[i] = 1;
    }
    prepare(n, &permanent[0], &base[0]);
    prepare_permanent(n, &erasures[0], &permanent[0], &base[0], &with_base[0]);
    Require(EquivalentLogs(with_base, expected),
        "permanent locator contribution mismatch");

    std::vector<uint8_t> dynamic_only(erasures);
    for (unsigned i = 0; i < n; ++i)
        if (permanent[i])
            dynamic_only[i] = 0;
    prepare_permanent(n, &dynamic_only[0], &permanent[0], &base[0],
        &with_dynamic_only[0]);
    Require(EquivalentLogs(with_dynamic_only, expected),
        "dynamic-only permanent locator contribution mismatch");
    return static_cast<uint64_t>(n) * (checked_direct ? 6 : 5) + sampled;
}

static unsigned AddModulo(unsigned a, unsigned b, unsigned modulus)
{
    return (a + b) % modulus;
}

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
                values[base + distance + i] = (a + modulus - b) % modulus;
            }
        }
    }
}

static std::vector<unsigned> DirectLocator(
    const std::vector<unsigned>& logarithms,
    const std::vector<uint8_t>& erasures,
    unsigned modulus)
{
    const unsigned n = static_cast<unsigned>(erasures.size());
    std::vector<unsigned> locator(n, 0);
    for (unsigned erased = 0; erased < n; ++erased)
    {
        if (!erasures[erased])
            continue;
        for (unsigned coordinate = 0; coordinate < n; ++coordinate)
        {
            locator[coordinate] = AddModulo(locator[coordinate],
                logarithms[coordinate ^ erased], modulus);
        }
    }
    return locator;
}

struct PermanentCandidateResult
{
    std::vector<unsigned> expected;
    std::vector<unsigned> coordinate_cache;
    std::vector<unsigned> transform_cache;
};

static PermanentCandidateResult EvaluatePermanentCandidates(
    const std::vector<unsigned>& logarithms,
    const std::vector<uint8_t>& permanent,
    const std::vector<uint8_t>& dynamic,
    unsigned field_order,
    unsigned modulus)
{
    const unsigned n = static_cast<unsigned>(permanent.size());
    Require(dynamic.size() == permanent.size(),
        "permanent candidate mask size mismatch");
    Require(field_order % n == 0,
        "permanent candidate parent must divide field order");

    std::vector<uint8_t> union_mask(n, 0);
    for (unsigned i = 0; i < n; ++i)
    {
        Require(!(permanent[i] && dynamic[i]),
            "permanent and dynamic masks must be disjoint");
        union_mask[i] = static_cast<uint8_t>(permanent[i] || dynamic[i]);
    }

    PermanentCandidateResult result;
    result.expected = DirectLocator(logarithms, union_mask, modulus);
    const std::vector<unsigned> permanent_coordinate =
        DirectLocator(logarithms, permanent, modulus);
    const std::vector<unsigned> dynamic_coordinate =
        DirectLocator(logarithms, dynamic, modulus);
    result.coordinate_cache.resize(n);
    for (unsigned i = 0; i < n; ++i)
    {
        result.coordinate_cache[i] = AddModulo(
            permanent_coordinate[i], dynamic_coordinate[i], modulus);
    }

    // Cache the permanent contribution after the first Walsh transform and
    // fixed-kernel multiplication.  This is the strongest transform-domain
    // candidate: it fuses the cache addition into the dynamic pointwise pass.
    std::vector<unsigned> kernel(logarithms);
    Walsh(kernel, modulus);
    const unsigned inverse_n = field_order / n;
    for (unsigned i = 0; i < n; ++i)
        kernel[i] = (kernel[i] * inverse_n) % modulus;

    std::vector<unsigned> permanent_frequency(n), dynamic_frequency(n);
    for (unsigned i = 0; i < n; ++i)
    {
        permanent_frequency[i] = permanent[i] ? 1u : 0u;
        dynamic_frequency[i] = dynamic[i] ? 1u : 0u;
    }
    Walsh(permanent_frequency, modulus);
    Walsh(dynamic_frequency, modulus);
    for (unsigned i = 0; i < n; ++i)
    {
        permanent_frequency[i] =
            (permanent_frequency[i] * kernel[i]) % modulus;
        dynamic_frequency[i] = AddModulo(
            (dynamic_frequency[i] * kernel[i]) % modulus,
            permanent_frequency[i], modulus);
    }
    Walsh(dynamic_frequency, modulus);
    result.transform_cache.swap(dynamic_frequency);
    return result;
}

static uint64_t CheckCandidateEquality(
    const std::vector<unsigned>& logarithms,
    const std::vector<uint8_t>& permanent,
    const std::vector<uint8_t>& dynamic,
    unsigned field_order,
    unsigned modulus)
{
    const PermanentCandidateResult result = EvaluatePermanentCandidates(
        logarithms, permanent, dynamic, field_order, modulus);
    Require(result.coordinate_cache == result.expected,
        "coordinate-domain permanent cache mismatch");
    Require(result.transform_cache == result.expected,
        "transform-domain permanent cache mismatch");
    return static_cast<uint64_t>(result.expected.size()) * 2u;
}

static uint64_t CheckExhaustiveGF4PermanentDecomposition()
{
    static const unsigned order = 16;
    static const unsigned modulus = order - 1;
    static const unsigned polynomial = 0x13;
    std::vector<unsigned> logarithms(order, 0);
    unsigned state = 1;
    for (unsigned exponent = 0; exponent < modulus; ++exponent)
    {
        logarithms[state] = exponent;
        state <<= 1;
        if (state >= order)
            state ^= polynomial;
    }

    uint64_t compared = 0;
    for (unsigned n = 2; n <= 8; n <<= 1)
    {
        const std::vector<unsigned> active_logs(logarithms.begin(),
            logarithms.begin() + n);
        uint32_t assignment_count = 1;
        for (unsigned i = 0; i < n; ++i)
            assignment_count *= 3;
        for (uint32_t assignment = 0; assignment < assignment_count;
             ++assignment)
        {
            uint32_t encoded = assignment;
            std::vector<uint8_t> permanent(n, 0), dynamic(n, 0);
            for (unsigned i = 0; i < n; ++i)
            {
                const unsigned role = encoded % 3u;
                encoded /= 3u;
                permanent[i] = static_cast<uint8_t>(role == 1);
                dynamic[i] = static_cast<uint8_t>(role == 2);
            }
            compared += CheckCandidateEquality(active_logs, permanent, dynamic,
                order, modulus);
        }
    }

    // N=16 already has every union subset checked by the normalization test.
    // Exercise deterministic decompositions here without turning this focused
    // suite into all 3^16 permanent/dynamic assignments.
    for (uint32_t seed = 0; seed < 4096; ++seed)
    {
        std::vector<uint8_t> permanent(order, 0), dynamic(order, 0);
        for (unsigned i = 0; i < order; ++i)
        {
            const unsigned role = Mix(seed * 31u + i * 0x9e37u) % 3u;
            permanent[i] = static_cast<uint8_t>(role == 1);
            dynamic[i] = static_cast<uint8_t>(role == 2);
        }
        compared += CheckCandidateEquality(logarithms, permanent, dynamic,
            order, modulus);
    }
    return compared;
}

struct StructuralCounts
{
    uint64_t cases;
    uint64_t sparse_additions_saved;

    StructuralCounts() : cases(0), sparse_additions_saved(0) {}
};

static StructuralCounts CheckPermanentCacheOperationCounts()
{
    StructuralCounts counts;
    for (unsigned n = 2; n <= 65536; n <<= 1)
    {
        unsigned log_n = 0;
        for (unsigned value = n; value > 1; value >>= 1)
            ++log_n;

        // One FWHT has N*log2(N) modular add/sub results.  Both dense cache
        // candidates retain both transforms and N pointwise multiplications,
        // then add one N-entry modular pass that the union-mask path avoids.
        const uint64_t union_modular = 2ull * n * log_n;
        const uint64_t coordinate_modular = union_modular + n;
        const uint64_t transform_modular = union_modular + n;
        Require(coordinate_modular > union_modular,
            "coordinate dense cache was not structurally dominated");
        Require(transform_modular > union_modular,
            "transform dense cache was not structurally dominated");

        const uint64_t permanent_count = std::max(1u, n / 8u);
        const uint64_t dynamic_count = std::min(7u, n - 1u);
        const uint64_t uncached_direct =
            (permanent_count + dynamic_count) * (n - 1u);
        const uint64_t cached_direct = dynamic_count * (n - 1u);
        Require(cached_direct < uncached_direct,
            "sparse permanent cache did not save direct additions");
        counts.sparse_additions_saved += uncached_direct - cached_direct;
        ++counts.cases;
    }
    return counts;
}

static uint64_t CheckSelectedErasureInvariant()
{
    uint64_t patterns = 0;
    for (unsigned k = 1; k <= 6; ++k)
    {
        for (unsigned r = 1; r <= 6; ++r)
        {
            const unsigned transmitted = k + r;
            const uint32_t pattern_count = 1u << transmitted;
            for (uint32_t missing = 0; missing < pattern_count; ++missing)
            {
                if (Popcount32(missing) > r)
                    continue;

                unsigned dynamic_erasure_count = 0;
                unsigned public_survivors_needed = k;
                for (unsigned original = 0; original < k; ++original)
                {
                    if ((missing >> original) & 1u)
                        ++dynamic_erasure_count;
                    else
                        --public_survivors_needed;
                }
                for (unsigned recovery = 0; recovery < r; ++recovery)
                {
                    const bool present =
                        ((missing >> (k + recovery)) & 1u) == 0;
                    if (!present)
                    {
                        ++dynamic_erasure_count;
                        continue;
                    }
                    if (public_survivors_needed != 0)
                        --public_survivors_needed;
                    else
                        ++dynamic_erasure_count; // deterministic virtual erasure
                }
                Require(public_survivors_needed == 0,
                    "valid pattern did not retain K public survivors");
                Require(dynamic_erasure_count == r,
                    "selected dynamic plus virtual erasures did not equal R");
                ++patterns;
            }
        }
    }
    return patterns;
}

static uint64_t CheckExhaustiveGF4Normalization()
{
    static const unsigned order = 16;
    static const unsigned modulus = order - 1;
    static const unsigned polynomial = 0x13;
    unsigned logarithm[order] = {};
    unsigned state = 1;
    for (unsigned exponent = 0; exponent < modulus; ++exponent)
    {
        logarithm[state] = exponent;
        state <<= 1;
        if (state >= order)
            state ^= polynomial;
    }

    uint64_t compared = 0;
    for (unsigned n = 2; n <= order; n <<= 1)
    {
        std::vector<unsigned> transformed_kernel(n);
        for (unsigned i = 1; i < n; ++i)
            transformed_kernel[i] = logarithm[i];
        Walsh(transformed_kernel, modulus);
        const unsigned inverse_n = order / n;
        for (unsigned i = 0; i < n; ++i)
            transformed_kernel[i] =
                (transformed_kernel[i] * inverse_n) % modulus;

        const uint32_t pattern_count = 1u << n;
        for (uint32_t pattern = 0; pattern < pattern_count; ++pattern)
        {
            std::vector<unsigned> actual(n);
            for (unsigned i = 0; i < n; ++i)
                actual[i] = (pattern >> i) & 1u;
            Walsh(actual, modulus);
            for (unsigned i = 0; i < n; ++i)
                actual[i] =
                    (actual[i] * transformed_kernel[i]) % modulus;
            Walsh(actual, modulus);

            for (unsigned coordinate = 0; coordinate < n; ++coordinate)
            {
                unsigned expected = 0;
                for (unsigned erased = 0; erased < n; ++erased)
                {
                    if (((pattern >> erased) & 1u) && erased != coordinate)
                    {
                        expected = AddModulo(
                            expected, logarithm[coordinate ^ erased], modulus);
                    }
                }
                Require(actual[coordinate] == expected,
                    "GF4 active Walsh normalization mismatch");
                ++compared;
            }
        }
    }
    return compared;
}

static uint64_t CheckExhaustiveGF8Active()
{
    uint64_t compared = 0;
    for (unsigned n = 2; n <= 16; n <<= 1)
    {
        const uint32_t pattern_count = 1u << n;
        std::vector<uint8_t> erasures(n), actual(n), expected(n);
        for (uint32_t pattern = 0; pattern < pattern_count; ++pattern)
        {
            for (unsigned i = 0; i < n; ++i)
                erasures[i] = static_cast<uint8_t>((pattern >> i) & 1u);
            leopard::ff8::PrepareDecodeWalshActive(
                n, &erasures[0], &actual[0]);

            for (unsigned coordinate = 0; coordinate < n; ++coordinate)
            {
                unsigned sum = 0;
                for (unsigned erased = 0; erased < n; ++erased)
                {
                    if (erasures[erased] && erased != coordinate)
                    {
                        sum = AddModulo(sum,
                            leopard::ff8::ElementLog(
                                static_cast<uint8_t>(coordinate ^ erased)),
                            leopard::ff8::kModulus);
                    }
                }
                expected[coordinate] = static_cast<uint8_t>(sum);
            }
            Require(EquivalentLogs(actual, expected),
                "GF8 exhaustive active/direct locator mismatch");
            compared += n;
        }
    }
    return compared;
}

template<class Ffe, class Prepare, class PreparePermanent, class Direct,
    class Active>
static uint64_t CheckPermanentTransitions(
    unsigned n,
    unsigned permanent_count,
    uint32_t seed,
    Prepare prepare,
    PreparePermanent prepare_permanent,
    Direct direct,
    Active active,
    Ffe (*element_log)(Ffe),
    unsigned field_order,
    unsigned modulus,
    bool include_dense)
{
    Require(permanent_count > 0 && permanent_count < n,
        "invalid permanent transition count");
    std::vector<unsigned> logarithms(n, 0);
    for (unsigned i = 1; i < n; ++i)
    {
        logarithms[i] = CanonicalLog(
            element_log(static_cast<Ffe>(i)));
    }

    std::vector<unsigned> order(n);
    for (unsigned i = 0; i < n; ++i)
        order[i] = i;
    std::sort(order.begin(), order.end(), [seed](unsigned a, unsigned b) {
        const uint32_t ma = Mix(a ^ seed);
        const uint32_t mb = Mix(b ^ seed);
        return ma == mb ? a < b : ma < mb;
    });

    std::vector<uint8_t> permanent(n, 0);
    for (unsigned i = 0; i < permanent_count; ++i)
        permanent[order[i]] = 1;

    std::vector<unsigned> available;
    for (size_t i = permanent_count; i < order.size(); ++i)
        available.push_back(order[i]);
    Require(available.size() >= 3, "insufficient dynamic transition space");

    std::vector<std::vector<uint8_t> > transitions;
    transitions.push_back(std::vector<uint8_t>(n, 0));

    std::vector<uint8_t> dynamic(n, 0);
    dynamic[available[0]] = 1; // add
    transitions.push_back(dynamic);
    const unsigned sparse_count = std::min<unsigned>(7, available.size());
    Require(available.size() > sparse_count,
        "insufficient space for a sparse swap transition");
    for (unsigned i = 1; i < sparse_count; ++i)
        dynamic[available[i]] = 1;
    transitions.push_back(dynamic); // multi-add
    dynamic[available[0]] = 0;
    transitions.push_back(dynamic); // remove
    dynamic[available[1]] = 0;
    dynamic[available[sparse_count]] = 1;
    transitions.push_back(dynamic); // swap
    if (include_dense)
    {
        dynamic.assign(n, 0);
        const unsigned dense_count = std::min<unsigned>(
            available.size(), std::max(8u, n / 4u));
        for (unsigned i = 0; i < dense_count; ++i)
            dynamic[available[i]] = 1;
        transitions.push_back(dynamic);
    }

    std::vector<Ffe> permanent_base(n);
    prepare(n, &permanent[0], &permanent_base[0]);
    uint64_t compared = 0;
    for (size_t transition = 0; transition < transitions.size(); ++transition)
    {
        const std::vector<uint8_t>& current = transitions[transition];
        const PermanentCandidateResult candidates =
            EvaluatePermanentCandidates(logarithms, permanent, current,
                field_order, modulus);
        Require(candidates.coordinate_cache == candidates.expected,
            "production transition coordinate candidate mismatch");
        Require(candidates.transform_cache == candidates.expected,
            "production transition transform candidate mismatch");

        std::vector<uint8_t> union_mask(permanent);
        for (unsigned i = 0; i < n; ++i)
            union_mask[i] = static_cast<uint8_t>(union_mask[i] || current[i]);
        std::vector<Ffe> active_actual(n), direct_actual(n), dynamic_actual(n),
            union_actual(n);
        active(n, &union_mask[0], &active_actual[0]);
        prepare_permanent(n, &current[0], &permanent[0], &permanent_base[0],
            &dynamic_actual[0]);
        prepare_permanent(n, &union_mask[0], &permanent[0],
            &permanent_base[0], &union_actual[0]);

        const unsigned total_count = static_cast<unsigned>(std::count(
            union_mask.begin(), union_mask.end(), uint8_t(1)));
        const bool check_full_direct =
            static_cast<uint64_t>(n) * total_count <= 4u * 1024u * 1024u;
        if (check_full_direct)
            direct(n, &union_mask[0], &direct_actual[0]);

        for (unsigned i = 0; i < n; ++i)
        {
            const unsigned expected = candidates.expected[i];
            Require(CanonicalLog(active_actual[i]) == expected,
                "production active transition mismatch");
            Require(CanonicalLog(dynamic_actual[i]) == expected,
                "production dynamic-only permanent transition mismatch");
            Require(CanonicalLog(union_actual[i]) == expected,
                "production union permanent transition mismatch");
            if (check_full_direct)
            {
                Require(CanonicalLog(direct_actual[i]) == expected,
                    "production direct transition mismatch");
            }
        }
        compared += static_cast<uint64_t>(n) *
            (check_full_direct ? 7u : 6u);
    }
    return compared;
}

static uint64_t CheckFF8()
{
    uint64_t compared = 0;
    for (unsigned n = 2; n <= leopard::ff8::kOrder; n <<= 1)
    {
        const unsigned counts[] = {
            0, 1, std::min(4u, n), std::min(5u, n), n / 8, n / 3, n / 2
        };
        for (unsigned i = 0; i < sizeof(counts) / sizeof(counts[0]); ++i)
            compared += CheckPattern<leopard::ff8::ffe_t>(
                n, counts[i], 0x81c5u + n * 17u + i,
                leopard::ff8::PrepareDecode,
                leopard::ff8::PrepareDecodeWithPermanent,
                leopard::ff8::PrepareDecodeDirect,
                leopard::ff8::PrepareDecodeWalshActive,
                leopard::ff8::PrepareDecodeWalshReference,
                leopard::ff8::ElementLog,
                leopard::ff8::kModulus);

        unsigned cutoff;
        if (n <= 8)
            cutoff = n;
        else if (n == 16)
            cutoff = 8;
        else if (n == 32 || n == 128)
            cutoff = 9;
        else if (n == 64)
            cutoff = 8;
        else
            cutoff = 7;
        Require(leopard::ff8::IsDirectLocatorPreferred(n, cutoff),
            "GF8 direct locator cutoff rejected");
        if (cutoff < n)
        {
            Require(!leopard::ff8::IsDirectLocatorPreferred(n, cutoff + 1),
                "GF8 dense locator cutoff accepted");
        }
    }

    const unsigned transition_sizes[] = { 16, 64, 256 };
    for (unsigned i = 0;
         i < sizeof(transition_sizes) / sizeof(transition_sizes[0]); ++i)
    {
        const unsigned n = transition_sizes[i];
        compared += CheckPermanentTransitions<leopard::ff8::ffe_t>(
            n, std::max(1u, n / 8u), 0x8ca11u + n,
            leopard::ff8::PrepareDecode,
            leopard::ff8::PrepareDecodeWithPermanent,
            leopard::ff8::PrepareDecodeDirect,
            leopard::ff8::PrepareDecodeWalshActive,
            leopard::ff8::ElementLog,
            leopard::ff8::kOrder,
            leopard::ff8::kModulus,
            true);
    }
    return compared;
}

static uint64_t CheckFF16()
{
    uint64_t compared = 0;
    for (unsigned n = 2; n <= leopard::ff16::kOrder; n <<= 1)
    {
        const unsigned counts[] = {
            0, 1, std::min(4u, n), std::min(5u, n), n / 16, n / 4, n / 2
        };
        for (unsigned i = 0; i < sizeof(counts) / sizeof(counts[0]); ++i)
            compared += CheckPattern<leopard::ff16::ffe_t>(
                n, counts[i], 0x16f00du + n * 31u + i,
                leopard::ff16::PrepareDecode,
                leopard::ff16::PrepareDecodeWithPermanent,
                leopard::ff16::PrepareDecodeDirect,
                leopard::ff16::PrepareDecodeWalshActive,
                leopard::ff16::PrepareDecodeWalshReference,
                leopard::ff16::ElementLog,
                leopard::ff16::kModulus);

        unsigned cutoff;
        if (n <= 32)
            cutoff = n;
        else if (n == 64)
            cutoff = 34;
        else if (n == 128)
            cutoff = 24;
        else if (n == 256 || n == 512)
            cutoff = 16;
        else
            cutoff = 14;
        Require(leopard::ff16::IsDirectLocatorPreferred(n, cutoff),
            "GF16 direct locator cutoff rejected");
        if (cutoff < n)
        {
            Require(!leopard::ff16::IsDirectLocatorPreferred(n, cutoff + 1),
                "GF16 dense locator cutoff accepted");
        }
    }

    const unsigned transition_sizes[] = { 32, 256, 2048, 65536 };
    for (unsigned i = 0;
         i < sizeof(transition_sizes) / sizeof(transition_sizes[0]); ++i)
    {
        const unsigned n = transition_sizes[i];
        // Keep the full-field direct oracle bounded while still exercising a
        // proper-parent dense permanent mask at the practical sizes.
        const unsigned permanent_count = n == 65536 ? 13u :
            std::max(1u, n / 8u);
        compared += CheckPermanentTransitions<leopard::ff16::ffe_t>(
            n, permanent_count, 0x16ca11u + n,
            leopard::ff16::PrepareDecode,
            leopard::ff16::PrepareDecodeWithPermanent,
            leopard::ff16::PrepareDecodeDirect,
            leopard::ff16::PrepareDecodeWalshActive,
            leopard::ff16::ElementLog,
            leopard::ff16::kOrder,
            leopard::ff16::kModulus,
            n <= 2048);
    }

    return compared;
}

static void CheckConcurrentSetup()
{
    const unsigned n = 1024;
    std::vector<uint8_t> erasures(n, 0), permanent(n, 0), dynamic(n, 0);
    for (unsigned i = 0; i < 19; ++i)
    {
        const unsigned coordinate = Mix(i + 91u) & (n - 1);
        erasures[coordinate] = 1;
        if ((i % 3u) == 0)
            permanent[coordinate] = 1;
        else
            dynamic[coordinate] = 1;
    }
    std::vector<leopard::ff16::ffe_t> expected(n), permanent_base(n);
    leopard::ff16::PrepareDecodeWalshReference(
        n, &erasures[0], &expected[0]);
    leopard::ff16::PrepareDecode(
        n, &permanent[0], &permanent_base[0]);

    std::atomic<unsigned> failures(0);
    std::vector<std::thread> threads;
    for (unsigned worker = 0; worker < 8; ++worker)
    {
        threads.push_back(std::thread([&]() {
            for (unsigned repeat = 0; repeat < 32; ++repeat)
            {
                std::vector<leopard::ff16::ffe_t> actual(n);
                if ((repeat & 1u) == 0)
                {
                    leopard::ff16::PrepareDecodeWithPermanent(
                        n, &dynamic[0], &permanent[0], &permanent_base[0],
                        &actual[0]);
                }
                else
                {
                    leopard::ff16::PrepareDecodeWithPermanent(
                        n, &erasures[0], &permanent[0], &permanent_base[0],
                        &actual[0]);
                }
                if (!EquivalentLogs(actual, expected))
                    ++failures;
            }
        }));
    }
    for (size_t i = 0; i < threads.size(); ++i)
        threads[i].join();
    Require(failures.load() == 0,
        "concurrent immutable permanent locator preparation mismatch");
}

} // namespace

int main()
{
    Require(leopard::ff8::Initialize(), "GF8 initialization");
    Require(leopard::ff16::Initialize(), "GF16 initialization");
    Require(CanonicalLog(std::numeric_limits<uint16_t>::max()) == 0,
        "GF16 logarithm sentinel canonicalization");

    const StructuralCounts structural = CheckPermanentCacheOperationCounts();
    const uint64_t selection_patterns = CheckSelectedErasureInvariant();
    const uint64_t compared = CheckExhaustiveGF4Normalization() +
        CheckExhaustiveGF4PermanentDecomposition() +
        CheckExhaustiveGF8Active() + CheckFF8() + CheckFF16();
    CheckConcurrentSetup();
    std::cout << "leopard2 locator tests passed: locator_entries_compared="
              << compared << " permanent_cache_structural_cases="
              << structural.cases << " sparse_additions_saved="
              << structural.sparse_additions_saved
              << " selection_patterns=" << selection_patterns
              << " concurrent_preparations=256" << std::endl;
    return 0;
}
