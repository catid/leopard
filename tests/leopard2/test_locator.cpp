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

template<class Ffe, class Prepare, class PreparePermanent, class Reference>
static uint64_t CheckPattern(
    unsigned n,
    unsigned erasure_count,
    uint32_t seed,
    Prepare prepare,
    PreparePermanent prepare_permanent,
    Reference reference)
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

    std::vector<Ffe> expected(n), actual(n), base(n), with_base(n);
    reference(n, &erasures[0], &expected[0]);
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
    return static_cast<uint64_t>(n) * 3;
}

static uint64_t CheckFF8()
{
    uint64_t compared = 0;
    for (unsigned n = 2; n <= leopard::ff8::kOrder; n <<= 1)
    {
        const unsigned counts[] = { 0, 1, n / 8, n / 3, n / 2 };
        for (unsigned i = 0; i < sizeof(counts) / sizeof(counts[0]); ++i)
            compared += CheckPattern<leopard::ff8::ffe_t>(
                n, counts[i], 0x81c5u + n * 17u + i,
                leopard::ff8::PrepareDecode,
                leopard::ff8::PrepareDecodeWithPermanent,
                leopard::ff8::PrepareDecodeWalshReference);
    }
    return compared;
}

static uint64_t CheckFF16()
{
    uint64_t compared = 0;
    for (unsigned n = 2; n <= 4096; n <<= 1)
    {
        const unsigned counts[] = { 0, 1, n / 16, n / 4, n / 2 };
        for (unsigned i = 0; i < sizeof(counts) / sizeof(counts[0]); ++i)
            compared += CheckPattern<leopard::ff16::ffe_t>(
                n, counts[i], 0x16f00du + n * 31u + i,
                leopard::ff16::PrepareDecode,
                leopard::ff16::PrepareDecodeWithPermanent,
                leopard::ff16::PrepareDecodeWalshReference);
    }

    // Exercise the dispatch crossover on a large active parent without making
    // the normal correctness test enumerate a quadratic direct pattern.
    compared += CheckPattern<leopard::ff16::ffe_t>(
        65536, 64, 0xf16f16u,
        leopard::ff16::PrepareDecode,
        leopard::ff16::PrepareDecodeWithPermanent,
        leopard::ff16::PrepareDecodeWalshReference);
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

    const uint64_t compared = CheckFF8() + CheckFF16();
    CheckConcurrentSetup();
    std::cout << "leopard2 locator tests passed: locator_entries_compared="
              << compared << " concurrent_preparations=256" << std::endl;
    return 0;
}
