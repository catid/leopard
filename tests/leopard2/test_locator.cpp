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

    return compared;
}

static void CheckConcurrentSetup()
{
    const unsigned n = 1024;
    std::vector<uint8_t> erasures(n, 0);
    for (unsigned i = 0; i < 19; ++i)
        erasures[Mix(i + 91u) & (n - 1)] = 1;
    std::vector<leopard::ff16::ffe_t> expected(n);
    leopard::ff16::PrepareDecodeWalshReference(
        n, &erasures[0], &expected[0]);

    std::atomic<unsigned> failures(0);
    std::vector<std::thread> threads;
    for (unsigned worker = 0; worker < 8; ++worker)
    {
        threads.push_back(std::thread([&]() {
            for (unsigned repeat = 0; repeat < 32; ++repeat)
            {
                std::vector<leopard::ff16::ffe_t> actual(n);
                leopard::ff16::PrepareDecode(n, &erasures[0], &actual[0]);
                if (!EquivalentLogs(actual, expected))
                    ++failures;
            }
        }));
    }
    for (size_t i = 0; i < threads.size(); ++i)
        threads[i].join();
    Require(failures.load() == 0, "concurrent locator preparation mismatch");
}

} // namespace

int main()
{
    Require(leopard::ff8::Initialize(), "GF8 initialization");
    Require(leopard::ff16::Initialize(), "GF16 initialization");
    Require(CanonicalLog(std::numeric_limits<uint16_t>::max()) == 0,
        "GF16 logarithm sentinel canonicalization");

    const uint64_t compared = CheckExhaustiveGF4Normalization() +
        CheckExhaustiveGF8Active() + CheckFF8() + CheckFF16();
    CheckConcurrentSetup();
    std::cout << "leopard2 locator tests passed: locator_entries_compared="
              << compared << " concurrent_preparations=256" << std::endl;
    return 0;
}
