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

// Default-off C8 experiment.  The root CMake graph does not reference this
// source, and no public profile or dispatcher path is changed by it.

#include "LeopardFF8.h"
#include "LeopardFF16.h"
#include "leopard2.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <new>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#if defined(__linux__)
#include <sched.h>
#endif

#if defined(_OPENMP)
#include <omp.h>
#endif

#ifndef LEO2_C8_SOURCE_SHA256
#define LEO2_C8_SOURCE_SHA256 "unbound"
#endif
#ifndef LEO2_C8_CORE_GIT_SHA
#define LEO2_C8_CORE_GIT_SHA "unbound"
#endif
#ifndef LEO2_C8_LIBRARY_SHA256
#define LEO2_C8_LIBRARY_SHA256 "unbound"
#endif
#ifndef LEO2_C8_SANITIZER_MODE
#define LEO2_C8_SANITIZER_MODE "none"
#endif

static bool C8TrackAllocations = false;
static uint64_t C8TrackedAllocations = 0;

#if !defined(LEO2_C8_DISABLE_GLOBAL_NEW_TRACKING)
#if defined(__GNUC__) || defined(__clang__)
#define C8_NOINLINE __attribute__((noinline))
#else
#define C8_NOINLINE
#endif
static C8_NOINLINE void* C8Allocate(size_t bytes)
{
    return malloc(bytes ? bytes : 1);
}
static C8_NOINLINE void C8Deallocate(void* pointer)
{
    free(pointer);
}
void* operator new(size_t bytes)
{
    if (C8TrackAllocations)
        ++C8TrackedAllocations;
    void* pointer = C8Allocate(bytes);
    if (!pointer)
        throw std::bad_alloc();
    return pointer;
}
void* operator new[](size_t bytes) { return ::operator new(bytes); }
void operator delete(void* pointer) noexcept { C8Deallocate(pointer); }
void operator delete[](void* pointer) noexcept { C8Deallocate(pointer); }
#define LEO2_C8_ALLOCATION_TRACKING "global-new"
#else
#define LEO2_C8_ALLOCATION_TRACKING "disabled-for-sanitizer"
#endif

namespace {

typedef std::chrono::steady_clock Clock;

static void Fail(const std::string& message)
{
    throw std::runtime_error(message);
}

static void Require(leo2_result result, const char* operation)
{
    if (result != LEO2_SUCCESS)
    {
        std::ostringstream stream;
        stream << operation << ": " << leo2_result_string(result);
        Fail(stream.str());
    }
}

static unsigned NextPow2(unsigned value)
{
    unsigned result = 1;
    while (result < value)
        result <<= 1;
    return result;
}

static uint64_t Fnv(uint64_t hash, const uint8_t* data, size_t bytes)
{
    for (size_t i = 0; i < bytes; ++i)
    {
        hash ^= data[i];
        hash *= UINT64_C(1099511628211);
    }
    return hash;
}

static const char* BackendName(leo2_backend backend)
{
    switch (backend)
    {
    case LEO2_BACKEND_SCALAR: return "scalar";
    case LEO2_BACKEND_SSSE3: return "ssse3";
    case LEO2_BACKEND_AVX2: return "avx2";
    case LEO2_BACKEND_NEON: return "neon";
    default: return "auto";
    }
}

template<class Element>
struct IndependentField;

template<>
struct IndependentField<uint8_t>
{
    uint8_t to_polynomial[256];
    uint8_t to_coordinate[256];

    IndependentField()
    {
        static const unsigned basis[8] = {
            1, 214, 152, 146, 86, 200, 88, 230
        };
        memset(to_coordinate, 0xff, sizeof(to_coordinate));
        for (unsigned coordinate = 0; coordinate < 256; ++coordinate)
        {
            unsigned polynomial = 0;
            for (unsigned bit = 0; bit < 8; ++bit)
                if (coordinate & (1U << bit))
                    polynomial ^= basis[bit];
            if (polynomial >= 256 || to_coordinate[polynomial] != 0xff)
                Fail("independent GF8 coordinate basis is not invertible");
            to_polynomial[coordinate] = static_cast<uint8_t>(polynomial);
            to_coordinate[polynomial] = static_cast<uint8_t>(coordinate);
        }
    }

    uint8_t Multiply(uint8_t left, uint8_t right) const
    {
        unsigned a = to_polynomial[left];
        unsigned b = to_polynomial[right];
        unsigned product = 0;
        while (b)
        {
            if (b & 1U) product ^= a;
            b >>= 1;
            a <<= 1;
        }
        for (int bit = 14; bit >= 8; --bit)
            if (product & (1U << bit))
                product ^= 0x11dU << (bit - 8);
        return to_coordinate[product];
    }

    uint8_t Inverse(uint8_t value) const
    {
        if (!value) Fail("independent GF8 inverse of zero");
        uint8_t result = 1;
        unsigned exponent = 254;
        while (exponent)
        {
            if (exponent & 1U) result = Multiply(result, value);
            exponent >>= 1;
            if (exponent) value = Multiply(value, value);
        }
        return result;
    }
};

template<>
struct IndependentField<uint16_t>
{
    uint16_t to_polynomial[65536];
    uint16_t to_coordinate[65536];

    IndependentField()
    {
        static const unsigned basis[16] = {
            0x0001, 0xACCA, 0x3C0E, 0x163E,
            0xC582, 0xED2E, 0x914C, 0x4012,
            0x6C98, 0x10D8, 0x6A72, 0xB900,
            0xFDB8, 0xFB34, 0xFF38, 0x991E
        };
        std::vector<uint8_t> seen(65536, 0);
        for (unsigned coordinate = 0; coordinate < 65536; ++coordinate)
        {
            unsigned polynomial = 0;
            for (unsigned bit = 0; bit < 16; ++bit)
                if (coordinate & (1U << bit)) polynomial ^= basis[bit];
            if (polynomial >= 65536 || seen[polynomial])
                Fail("independent GF16 coordinate basis is not invertible");
            to_polynomial[coordinate] = static_cast<uint16_t>(polynomial);
            to_coordinate[polynomial] = static_cast<uint16_t>(coordinate);
            seen[polynomial] = 1;
        }
    }

    uint16_t Multiply(uint16_t left, uint16_t right) const
    {
        uint64_t a = to_polynomial[left];
        uint64_t b = to_polynomial[right];
        uint64_t product = 0;
        while (b)
        {
            if (b & 1U) product ^= a;
            b >>= 1;
            a <<= 1;
        }
        for (int bit = 30; bit >= 16; --bit)
            if (product & (UINT64_C(1) << bit))
                product ^= UINT64_C(0x1002d) << (bit - 16);
        return to_coordinate[static_cast<uint16_t>(product)];
    }

    uint16_t Inverse(uint16_t value) const
    {
        if (!value) Fail("independent GF16 inverse of zero");
        uint16_t result = 1;
        unsigned exponent = 65534;
        while (exponent)
        {
            if (exponent & 1U) result = Multiply(result, value);
            exponent >>= 1;
            if (exponent) value = Multiply(value, value);
        }
        return result;
    }
};

template<class Element>
struct ProductionField;

template<>
struct ProductionField<uint8_t>
{
    static uint8_t Multiply(uint8_t a, uint8_t b)
        { return leopard::ff8::MultiplyElements(a, b); }
    static uint8_t Inverse(uint8_t value)
        { return leopard::ff8::InverseElement(value); }
    static uint8_t Log(uint8_t value)
        { return leopard::ff8::ElementLog(value); }
    static uint8_t FromLog(uint8_t value)
        { return leopard::ff8::MultiplyLogElement(1, value); }
    static void MultiplyBytes(void* output, const void* input,
                              uint8_t log, uint64_t bytes)
        { leopard::ff8::MultiplyBytes(output, input, log, bytes); }
    static void MultiplyAddBytes(void* output, const void* input,
                                 uint8_t log, uint64_t bytes)
        { leopard::ff8::MultiplyAddBytes(output, input, log, bytes); }
    static leo2_field Field() { return LEO2_FIELD_GF8; }
    static const char* Name() { return "gf8"; }
};

template<>
struct ProductionField<uint16_t>
{
    static uint16_t Multiply(uint16_t a, uint16_t b)
        { return leopard::ff16::MultiplyElements(a, b); }
    static uint16_t Inverse(uint16_t value)
        { return leopard::ff16::InverseElement(value); }
    static uint16_t Log(uint16_t value)
        { return leopard::ff16::ElementLog(value); }
    static uint16_t FromLog(uint16_t value)
        { return leopard::ff16::MultiplyLogElement(1, value); }
    static void MultiplyBytes(void* output, const void* input,
                              uint16_t log, uint64_t bytes)
        { leopard::ff16::MultiplyBytes(output, input, log, bytes); }
    static void MultiplyAddBytes(void* output, const void* input,
                                 uint16_t log, uint64_t bytes)
        { leopard::ff16::MultiplyAddBytes(output, input, log, bytes); }
    static leo2_field Field() { return LEO2_FIELD_GF16; }
    static const char* Name() { return "gf16"; }
};

template<class Element>
class ExactPlan
{
public:
    unsigned k;
    unsigned r;
    std::vector<Element> logs;

    ExactPlan(unsigned original_count, unsigned recovery_count)
        : k(original_count), r(recovery_count), logs(k * r)
    {
        if (!k || !r || static_cast<uint64_t>(k) + r >
                (sizeof(Element) == 1 ? 256U : 65536U))
            Fail("invalid C8 exact plan geometry");
        std::vector<Element> weights(k);
        for (unsigned i = 0; i < k; ++i)
        {
            const Element point = static_cast<Element>(r + i);
            Element denominator = 1;
            for (unsigned other = 0; other < k; ++other)
            {
                if (other == i) continue;
                denominator = ProductionField<Element>::Multiply(
                    denominator, static_cast<Element>(point ^ (r + other)));
            }
            weights[i] = ProductionField<Element>::Inverse(denominator);
        }
        for (unsigned parity = 0; parity < r; ++parity)
        {
            const Element output_point = static_cast<Element>(parity);
            Element vanishing = 1;
            for (unsigned i = 0; i < k; ++i)
                vanishing = ProductionField<Element>::Multiply(
                    vanishing, static_cast<Element>(output_point ^ (r + i)));
            for (unsigned i = 0; i < k; ++i)
            {
                Element coefficient = ProductionField<Element>::Multiply(
                    vanishing, ProductionField<Element>::Inverse(
                        static_cast<Element>(output_point ^ (r + i))));
                coefficient = ProductionField<Element>::Multiply(
                    coefficient, weights[i]);
                if (!coefficient)
                    Fail("C8 exact generator unexpectedly contains zero");
                logs[parity * k + i] = ProductionField<Element>::Log(coefficient);
            }
        }
    }

    Element Coefficient(unsigned parity, unsigned source) const
    {
        return ProductionField<Element>::FromLog(logs[parity * k + source]);
    }

    void Encode(uint64_t bytes, const void* const* source,
                void* const* parity) const
    {
        for (unsigned output_index = 0; output_index < r; ++output_index)
        {
            void* output = parity[output_index];
            const Element first = logs[output_index * k];
            if (first == 0)
                memcpy(output, source[0], static_cast<size_t>(bytes));
            else
                ProductionField<Element>::MultiplyBytes(
                    output, source[0], first, bytes);
            for (unsigned input_index = 1; input_index < k; ++input_index)
                ProductionField<Element>::MultiplyAddBytes(
                    output, source[input_index],
                    logs[output_index * k + input_index], bytes);
        }
    }
};

template<class Element>
static std::vector<Element> IndependentRows(
    const IndependentField<Element>& field, unsigned k, unsigned r)
{
    std::vector<Element> weights(k);
    for (unsigned i = 0; i < k; ++i)
    {
        const Element point = static_cast<Element>(r + i);
        Element denominator = 1;
        for (unsigned other = 0; other < k; ++other)
            if (other != i)
                denominator = field.Multiply(
                    denominator, static_cast<Element>(point ^ (r + other)));
        weights[i] = field.Inverse(denominator);
    }
    std::vector<Element> rows(k * r);
    for (unsigned parity = 0; parity < r; ++parity)
    {
        const Element output_point = static_cast<Element>(parity);
        Element vanishing = 1;
        for (unsigned i = 0; i < k; ++i)
            vanishing = field.Multiply(
                vanishing, static_cast<Element>(output_point ^ (r + i)));
        for (unsigned i = 0; i < k; ++i)
            rows[parity * k + i] = field.Multiply(
                field.Multiply(vanishing, field.Inverse(
                    static_cast<Element>(output_point ^ (r + i)))), weights[i]);
    }
    return rows;
}

template<class Element>
static Element LoadSymbol(const uint8_t* data, size_t bytes, size_t symbol)
{
    if (sizeof(Element) == 1)
        return static_cast<Element>(data[symbol]);
    const size_t tile_symbol = symbol % 32;
    const size_t tile_start = (symbol / 32) * 64;
    const size_t remaining = bytes - tile_start;
    const size_t symbols = std::min<size_t>(32, remaining / 2);
    return static_cast<Element>(data[tile_start + tile_symbol] |
        (static_cast<unsigned>(data[tile_start + symbols + tile_symbol]) << 8));
}

template<class Element>
static void StoreSymbol(uint8_t* data, size_t bytes, size_t symbol, Element value)
{
    if (sizeof(Element) == 1)
    {
        data[symbol] = static_cast<uint8_t>(value);
        return;
    }
    const size_t tile_symbol = symbol % 32;
    const size_t tile_start = (symbol / 32) * 64;
    const size_t remaining = bytes - tile_start;
    const size_t symbols = std::min<size_t>(32, remaining / 2);
    data[tile_start + tile_symbol] = static_cast<uint8_t>(value);
    data[tile_start + symbols + tile_symbol] = static_cast<uint8_t>(value >> 8);
}

template<class Element>
static void IndependentEncode(
    const IndependentField<Element>& field,
    const std::vector<Element>& rows, unsigned k, unsigned r,
    size_t bytes, const std::vector<std::vector<uint8_t> >& source,
    std::vector<std::vector<uint8_t> >& parity)
{
    const size_t symbols = bytes / sizeof(Element);
    for (unsigned output = 0; output < r; ++output)
        for (size_t symbol = 0; symbol < symbols; ++symbol)
        {
            Element value = 0;
            for (unsigned input = 0; input < k; ++input)
                value ^= field.Multiply(rows[output * k + input],
                                        LoadSymbol<Element>(&source[input][0], bytes, symbol));
            StoreSymbol<Element>(&parity[output][0], bytes, symbol, value);
        }
}

struct Correctness
{
    uint64_t cases;
    uint64_t coefficients;
    uint64_t byte_comparisons;
    uint64_t legacy_identity_cases;
    uint64_t legacy_difference_cases;
    uint64_t hot_path_allocations;
    uint64_t digest;
};

template<class Element>
static void RunCase(unsigned k, unsigned r, unsigned case_index,
                    const std::vector<size_t>& byte_counts,
                    leo2_context* context, Correctness& result)
{
    ExactPlan<Element> exact(k, r);
    IndependentField<Element> independent;
    const std::vector<Element> rows = IndependentRows(independent, k, r);
    for (unsigned parity = 0; parity < r; ++parity)
        for (unsigned source = 0; source < k; ++source)
        {
            if (exact.Coefficient(parity, source) != rows[parity * k + source])
                Fail("production/independent C8 coefficient mismatch");
            ++result.coefficients;
        }

    leo2_codec_options options = {};
    options.struct_size = sizeof(options);
    leo2_codec* legacy = NULL;
    Require(leo2_codec_create(context, k, r, LEO2_PROFILE_LEGACY_HIGH_V1,
                              ProductionField<Element>::Field(), &options, &legacy),
            "legacy codec create");
    const size_t legacy_alignment = leo2_scratch_alignment();
    bool saw_legacy_difference = false;

    for (size_t byte_index = 0; byte_index < byte_counts.size(); ++byte_index)
    {
        const size_t bytes = byte_counts[byte_index];
        if (bytes % sizeof(Element))
            Fail("GF16 correctness size is not symbol aligned");
        size_t legacy_scratch_size = 0;
        Require(leo2_encode_scratch_size(legacy, bytes, &legacy_scratch_size),
                "legacy scratch query");
        std::vector<std::vector<uint8_t> > source(k, std::vector<uint8_t>(bytes));
        std::vector<std::vector<uint8_t> > observed(r, std::vector<uint8_t>(bytes));
        std::vector<std::vector<uint8_t> > expected(r, std::vector<uint8_t>(bytes));
        std::vector<std::vector<uint8_t> > legacy_output(r, std::vector<uint8_t>(bytes));
        std::vector<const void*> source_ptr(k);
        std::vector<void*> observed_ptr(r), legacy_ptr(r);
        for (unsigned input = 0; input < k; ++input)
        {
            for (size_t byte = 0; byte < bytes; ++byte)
                source[input][byte] = static_cast<uint8_t>(
                    input * 97U + byte * 53U + case_index * 29U + byte_index + 1U);
            source_ptr[input] = &source[input][0];
        }
        for (unsigned output = 0; output < r; ++output)
        {
            observed_ptr[output] = &observed[output][0];
            legacy_ptr[output] = &legacy_output[output][0];
        }
        C8TrackedAllocations = 0;
        C8TrackAllocations = true;
        exact.Encode(bytes, &source_ptr[0], &observed_ptr[0]);
        C8TrackAllocations = false;
        result.hot_path_allocations += C8TrackedAllocations;
        IndependentEncode(independent, rows, k, r, bytes, source, expected);
        for (unsigned output = 0; output < r; ++output)
        {
            if (observed[output] != expected[output])
                Fail("C8 byte execution differs from independent field oracle");
            result.byte_comparisons += bytes;
            result.digest = Fnv(result.digest, &observed[output][0], bytes);
        }

        void* legacy_scratch = NULL;
        if (legacy_scratch_size && posix_memalign(&legacy_scratch, legacy_alignment,
                                                  legacy_scratch_size) != 0)
            throw std::bad_alloc();
        Require(leo2_encode(legacy, bytes, &source_ptr[0], &legacy_ptr[0],
                            legacy_scratch, legacy_scratch_size),
                "legacy encode");
        free(legacy_scratch);
        const bool identity = (r & (r - 1)) == 0 && k + r == NextPow2(k + r);
        bool equal = true;
        for (unsigned output = 0; output < r; ++output)
            equal = equal && observed[output] == legacy_output[output];
        if (identity && !equal)
            Fail("C8 full-parent executable legacy wire identity failed");
        if (!identity && !equal)
            saw_legacy_difference = true;
        if (identity) ++result.legacy_identity_cases;
    }
    if ((r & (r - 1)) != 0 || k + r != NextPow2(k + r))
    {
        if (!saw_legacy_difference)
            Fail("C8 nonidentity geometry produced no executable wire witness");
        ++result.legacy_difference_cases;
    }
    leo2_codec_destroy(legacy);
    ++result.cases;
}

static Correctness RunCorrectness(leo2_context* context)
{
    Correctness result = { 0, 0, 0, 0, 0, 0, UINT64_C(1469598103934665603) };
    static const unsigned gf8_cases[][2] = {
        {1,1}, {3,1}, {5,3}, {12,4}, {24,8}, {31,3}, {48,16},
        {60,4}, {63,1}, {96,32}, {120,7}, {120,8}, {120,9},
        {192,64}, {224,31}
    };
    static const unsigned gf16_cases[][2] = {
        {255,1}, {500,3}, {1000,17}, {1000,33}
    };
    const size_t gf8_bytes_array[] = { 1, 7, 31, 64, 65, 257 };
    const size_t gf16_bytes_array[] = { 2, 6, 62, 64, 66, 130 };
    const std::vector<size_t> gf8_bytes(
        gf8_bytes_array, gf8_bytes_array + sizeof(gf8_bytes_array) / sizeof(size_t));
    const std::vector<size_t> gf16_bytes(
        gf16_bytes_array, gf16_bytes_array + sizeof(gf16_bytes_array) / sizeof(size_t));
    for (size_t i = 0; i < sizeof(gf8_cases) / sizeof(gf8_cases[0]); ++i)
        RunCase<uint8_t>(gf8_cases[i][0], gf8_cases[i][1],
                         static_cast<unsigned>(i), gf8_bytes, context, result);
    for (size_t i = 0; i < sizeof(gf16_cases) / sizeof(gf16_cases[0]); ++i)
        RunCase<uint16_t>(gf16_cases[i][0], gf16_cases[i][1],
                          static_cast<unsigned>(i + 100), gf16_bytes, context, result);
    if (result.hot_path_allocations != 0)
        Fail("C8 exact execution allocated in the hot path");
    return result;
}

struct Summary
{
    double median_us;
    double mad_us;
};

static Summary Summarize(std::vector<double> values)
{
    if (values.empty()) Fail("empty timing sample");
    std::sort(values.begin(), values.end());
    const double median = values[values.size() / 2];
    std::vector<double> deviations;
    for (size_t i = 0; i < values.size(); ++i)
        deviations.push_back(std::fabs(values[i] - median));
    std::sort(deviations.begin(), deviations.end());
    Summary result = { median, deviations[deviations.size() / 2] };
    return result;
}

struct BenchmarkSpec
{
    bool gf16;
    unsigned k;
    unsigned r;
    size_t bytes;
    unsigned batch;
    unsigned reuse;
};

struct BenchmarkResult
{
    BenchmarkSpec spec;
    Summary exact;
    Summary padded;
    Summary exact_setup;
    Summary padded_setup;
    size_t exact_table_bytes;
    size_t padded_scratch_bytes;
    uint64_t exact_logical_bytes;
    uint64_t digest;
    std::vector<double> exact_samples;
    std::vector<double> padded_samples;
    std::vector<double> exact_setup_samples;
    std::vector<double> padded_setup_samples;
};

template<class Element>
static BenchmarkResult RunBenchmarkCell(leo2_context* context,
                                        const BenchmarkSpec& spec)
{
    std::vector<double> exact_setup_samples, padded_setup_samples;
    for (unsigned sample = 0; sample < 7; ++sample)
    {
        Clock::time_point begin = Clock::now();
        ExactPlan<Element> plan(spec.k, spec.r);
        Clock::time_point middle = Clock::now();
        leo2_codec_options options = {};
        options.struct_size = sizeof(options);
        leo2_codec* codec = NULL;
        Require(leo2_codec_create(context, spec.k, spec.r,
                                  LEO2_PROFILE_LEGACY_HIGH_V1,
                                  ProductionField<Element>::Field(),
                                  &options, &codec),
                "benchmark codec setup");
        Clock::time_point end = Clock::now();
        exact_setup_samples.push_back(
            std::chrono::duration<double, std::micro>(middle - begin).count());
        padded_setup_samples.push_back(
            std::chrono::duration<double, std::micro>(end - middle).count());
        leo2_codec_destroy(codec);
    }

    ExactPlan<Element> exact(spec.k, spec.r);
    leo2_codec_options options = {};
    options.struct_size = sizeof(options);
    leo2_codec* codec = NULL;
    Require(leo2_codec_create(context, spec.k, spec.r,
                              LEO2_PROFILE_LEGACY_HIGH_V1,
                              ProductionField<Element>::Field(), &options, &codec),
            "benchmark codec create");
    size_t scratch_size = 0;
    Require(leo2_encode_scratch_size(codec, spec.bytes, &scratch_size),
            "benchmark scratch query");
    const size_t scratch_alignment = leo2_scratch_alignment();
    void* scratch = NULL;
    if (scratch_size && posix_memalign(&scratch, scratch_alignment, scratch_size) != 0)
        throw std::bad_alloc();

    const size_t stripes = spec.batch;
    std::vector<std::vector<std::vector<uint8_t> > > source(
        stripes, std::vector<std::vector<uint8_t> >(
            spec.k, std::vector<uint8_t>(spec.bytes)));
    std::vector<std::vector<std::vector<uint8_t> > > exact_output(
        stripes, std::vector<std::vector<uint8_t> >(
            spec.r, std::vector<uint8_t>(spec.bytes)));
    std::vector<std::vector<std::vector<uint8_t> > > padded_output = exact_output;
    std::vector<std::vector<const void*> > source_ptr(
        stripes, std::vector<const void*>(spec.k));
    std::vector<std::vector<void*> > exact_ptr(
        stripes, std::vector<void*>(spec.r));
    std::vector<std::vector<void*> > padded_ptr = exact_ptr;
    for (size_t stripe = 0; stripe < stripes; ++stripe)
    {
        for (unsigned input = 0; input < spec.k; ++input)
        {
            for (size_t byte = 0; byte < spec.bytes; ++byte)
                source[stripe][input][byte] = static_cast<uint8_t>(
                    stripe * 31U + input * 17U + byte * 13U + 1U);
            source_ptr[stripe][input] = &source[stripe][input][0];
        }
        for (unsigned output = 0; output < spec.r; ++output)
        {
            exact_ptr[stripe][output] = &exact_output[stripe][output][0];
            padded_ptr[stripe][output] = &padded_output[stripe][output][0];
        }
    }
    const unsigned warmups = 2;
    for (unsigned warmup = 0; warmup < warmups; ++warmup)
        for (size_t stripe = 0; stripe < stripes; ++stripe)
        {
            exact.Encode(spec.bytes, &source_ptr[stripe][0], &exact_ptr[stripe][0]);
            Require(leo2_encode(codec, spec.bytes, &source_ptr[stripe][0],
                                &padded_ptr[stripe][0], scratch, scratch_size),
                    "benchmark padded warmup");
        }

    std::vector<double> exact_samples, padded_samples;
    for (unsigned sample = 0; sample < 11; ++sample)
    {
        double elapsed[2] = { 0, 0 };
        for (unsigned order = 0; order < 2; ++order)
        {
            const bool run_exact = ((sample + order) & 1U) == 0;
            Clock::time_point begin = Clock::now();
            for (unsigned reuse = 0; reuse < spec.reuse; ++reuse)
                for (size_t stripe = 0; stripe < stripes; ++stripe)
                {
                    if (run_exact)
                        exact.Encode(spec.bytes, &source_ptr[stripe][0],
                                     &exact_ptr[stripe][0]);
                    else
                        Require(leo2_encode(
                            codec, spec.bytes, &source_ptr[stripe][0],
                            &padded_ptr[stripe][0], scratch, scratch_size),
                            "benchmark padded encode");
                }
            Clock::time_point end = Clock::now();
            elapsed[run_exact ? 0 : 1] =
                std::chrono::duration<double, std::micro>(end - begin).count() /
                (spec.reuse * stripes);
        }
        exact_samples.push_back(elapsed[0]);
        padded_samples.push_back(elapsed[1]);
    }
    uint64_t digest = UINT64_C(1469598103934665603);
    for (size_t stripe = 0; stripe < stripes; ++stripe)
        for (unsigned output = 0; output < spec.r; ++output)
        {
            digest = Fnv(digest, &exact_output[stripe][output][0], spec.bytes);
            digest = Fnv(digest, &padded_output[stripe][output][0], spec.bytes);
        }
    BenchmarkResult result = {
        spec, Summarize(exact_samples), Summarize(padded_samples),
        Summarize(exact_setup_samples), Summarize(padded_setup_samples),
        spec.k * spec.r * sizeof(Element), scratch_size,
        static_cast<uint64_t>(spec.r) *
            (1U + 3U * (spec.k - 1U)) * spec.bytes,
        digest, exact_samples, padded_samples,
        exact_setup_samples, padded_setup_samples
    };
    free(scratch);
    leo2_codec_destroy(codec);
    return result;
}

static std::vector<BenchmarkSpec> BenchmarkSpecs()
{
    const BenchmarkSpec values[] = {
        {false,31,2,64,8,8}, {false,31,2,1024,8,4},
        {false,31,3,64,8,8}, {false,31,3,1024,8,4},
        {false,31,3,65536,1,1},
        {false,31,4,64,8,8}, {false,31,4,1024,8,4},
        {false,60,4,1024,8,4}, {false,60,4,65536,1,1},
        {false,120,7,64,8,8},
        {false,120,7,1024,8,4}, {false,120,7,65536,1,1},
        {false,120,8,1024,8,4}, {false,120,9,1024,8,4},
        {false,224,31,64,8,4}, {false,224,31,1024,1,1},
        {true,500,2,64,8,4}, {true,500,2,1024,1,1},
        {true,500,3,64,8,4}, {true,500,3,1024,1,1},
        {true,500,4,64,8,4}, {true,500,4,1024,1,1},
        {true,500,5,64,8,4}, {true,500,5,1024,1,1},
        {true,1000,17,64,1,1}, {true,1000,17,1024,1,1},
        {true,1000,33,64,1,1}
    };
    return std::vector<BenchmarkSpec>(
        values, values + sizeof(values) / sizeof(values[0]));
}

static std::vector<int> Affinity()
{
    std::vector<int> result;
#if defined(__linux__)
    cpu_set_t set;
    CPU_ZERO(&set);
    if (sched_getaffinity(0, sizeof(set), &set) == 0)
        for (int cpu = 0; cpu < CPU_SETSIZE; ++cpu)
            if (CPU_ISSET(cpu, &set)) result.push_back(cpu);
#endif
    return result;
}

static void WriteJson(
    const std::string& path, const std::string& mode, const std::string& label,
    leo2_context* context, const Correctness* correctness,
    const std::vector<BenchmarkResult>& benchmarks)
{
    std::ofstream file;
    std::ostream* output = &std::cout;
    if (path != "-")
    {
        file.open(path.c_str(), std::ios::binary | std::ios::trunc);
        if (!file) Fail("cannot open output JSON");
        output = &file;
    }
    const std::vector<int> affinity = Affinity();
    *output << std::fixed << std::setprecision(6)
        << "{\"schema\":\"leopard2-c8-executable-v1\""
        << ",\"mode\":\"" << mode << "\""
        << ",\"backend_label\":\"" << label << "\""
        << ",\"runtime_backend\":\"" << BackendName(leo2_context_backend(context)) << "\""
        << ",\"source_sha256\":\"" << LEO2_C8_SOURCE_SHA256 << "\""
        << ",\"core_git_sha\":\"" << LEO2_C8_CORE_GIT_SHA << "\""
        << ",\"library_sha256\":\"" << LEO2_C8_LIBRARY_SHA256 << "\""
        << ",\"sanitizer\":\"" << LEO2_C8_SANITIZER_MODE << "\""
        << ",\"allocation_tracking\":\"" << LEO2_C8_ALLOCATION_TRACKING << "\""
        << ",\"profile\":\"exact_high_prefix_v1_candidate\""
        << ",\"default_off\":true"
        << ",\"openmp_max_threads\":";
#if defined(_OPENMP)
    *output << omp_get_max_threads();
#else
    *output << 1;
#endif
    *output << ",\"affinity\":[";
    for (size_t i = 0; i < affinity.size(); ++i)
    {
        if (i) *output << ',';
        *output << affinity[i];
    }
    *output << ']';
    if (correctness)
    {
        *output << ",\"correctness\":{"
            << "\"cases\":" << correctness->cases
            << ",\"coefficients\":" << correctness->coefficients
            << ",\"byte_comparisons\":" << correctness->byte_comparisons
            << ",\"legacy_identity_cases\":" << correctness->legacy_identity_cases
            << ",\"legacy_difference_cases\":" << correctness->legacy_difference_cases
            << ",\"hot_path_allocations\":" << correctness->hot_path_allocations
            << ",\"digest\":\"0x" << std::hex << correctness->digest << std::dec << "\"}";
    }
    *output << ",\"benchmarks\":[";
    for (size_t i = 0; i < benchmarks.size(); ++i)
    {
        if (i) *output << ',';
        const BenchmarkResult& row = benchmarks[i];
        const double ratio = row.padded.median_us / row.exact.median_us;
        const double credible = (row.padded.median_us - row.padded.mad_us) /
            (row.exact.median_us + row.exact.mad_us) - 1.;
        *output << "{\"field\":\"" << (row.spec.gf16 ? "gf16" : "gf8") << "\""
            << ",\"K\":" << row.spec.k << ",\"R\":" << row.spec.r
            << ",\"bytes\":" << row.spec.bytes
            << ",\"batch\":" << row.spec.batch << ",\"reuse\":" << row.spec.reuse
            << ",\"exact_median_us\":" << row.exact.median_us
            << ",\"exact_mad_us\":" << row.exact.mad_us
            << ",\"padded_median_us\":" << row.padded.median_us
            << ",\"padded_mad_us\":" << row.padded.mad_us
            << ",\"ratio\":" << ratio
            << ",\"credible_gain\":" << credible
            << ",\"exact_setup_median_us\":" << row.exact_setup.median_us
            << ",\"exact_setup_mad_us\":" << row.exact_setup.mad_us
            << ",\"padded_setup_median_us\":" << row.padded_setup.median_us
            << ",\"padded_setup_mad_us\":" << row.padded_setup.mad_us
            << ",\"exact_table_bytes\":" << row.exact_table_bytes
            << ",\"exact_scratch_bytes\":0"
            << ",\"padded_scratch_bytes\":" << row.padded_scratch_bytes
            << ",\"exact_logical_bytes\":" << row.exact_logical_bytes
            << ",\"digest\":\"0x" << std::hex << row.digest << std::dec << "\""
            << ",\"exact_samples_us\":[";
        for (size_t sample = 0; sample < row.exact_samples.size(); ++sample)
        {
            if (sample) *output << ',';
            *output << row.exact_samples[sample];
        }
        *output << "],\"padded_samples_us\":[";
        for (size_t sample = 0; sample < row.padded_samples.size(); ++sample)
        {
            if (sample) *output << ',';
            *output << row.padded_samples[sample];
        }
        *output << "],\"exact_setup_samples_us\":[";
        for (size_t sample = 0; sample < row.exact_setup_samples.size(); ++sample)
        {
            if (sample) *output << ',';
            *output << row.exact_setup_samples[sample];
        }
        *output << "],\"padded_setup_samples_us\":[";
        for (size_t sample = 0; sample < row.padded_setup_samples.size(); ++sample)
        {
            if (sample) *output << ',';
            *output << row.padded_setup_samples[sample];
        }
        *output << "]}";
    }
    *output << "]}\n";
}

struct Options
{
    std::string mode;
    std::string backend_label;
    std::string output;
    Options() : mode("correctness"), backend_label("auto"), output("-") {}
};

static Options Parse(int argc, char** argv)
{
    Options options;
    for (int i = 1; i < argc; ++i)
    {
        const std::string argument(argv[i]);
        if ((argument == "--mode" || argument == "--backend-label" ||
             argument == "--output") && i + 1 >= argc)
            Fail("missing option value");
        if (argument == "--mode") options.mode = argv[++i];
        else if (argument == "--backend-label") options.backend_label = argv[++i];
        else if (argument == "--output") options.output = argv[++i];
        else if (argument == "--help")
        {
            std::cout << "Usage: c8_exact_high --mode correctness|benchmark|all "
                         "--backend-label NAME --output FILE\n";
            exit(0);
        }
        else Fail("unknown option: " + argument);
    }
    if (options.mode != "correctness" && options.mode != "benchmark" &&
        options.mode != "all")
        Fail("invalid --mode");
    return options;
}

} // namespace

int main(int argc, char** argv)
{
    try
    {
        const Options options = Parse(argc, argv);
        leo2_context_options context_options = {};
        context_options.struct_size = sizeof(context_options);
        context_options.backend = LEO2_BACKEND_AUTO;
        context_options.thread_count = 1;
        leo2_context* context = NULL;
        Require(leo2_context_create(&context_options, &context), "context create");
        Correctness correctness;
        Correctness* correctness_pointer = NULL;
        if (options.mode == "correctness" || options.mode == "all")
        {
            correctness = RunCorrectness(context);
            correctness_pointer = &correctness;
        }
        std::vector<BenchmarkResult> benchmark_results;
        if (options.mode == "benchmark" || options.mode == "all")
        {
            const std::vector<BenchmarkSpec> specs = BenchmarkSpecs();
            for (size_t i = 0; i < specs.size(); ++i)
            {
                if (specs[i].gf16)
                    benchmark_results.push_back(
                        RunBenchmarkCell<uint16_t>(context, specs[i]));
                else
                    benchmark_results.push_back(
                        RunBenchmarkCell<uint8_t>(context, specs[i]));
            }
        }
        WriteJson(options.output, options.mode, options.backend_label, context,
                  correctness_pointer, benchmark_results);
        leo2_context_destroy(context);
        return 0;
    }
    catch (const std::exception& error)
    {
        std::cerr << "C8 failure: " << error.what() << '\n';
        return 1;
    }
}
