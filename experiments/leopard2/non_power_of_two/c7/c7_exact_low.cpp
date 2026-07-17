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

// Default-off C7 exact-low experiment.  The root CMake graph does not compile
// this file and production constructors continue to reject profile enum 3.

#include "LeopardFF8.h"
#include "LeopardFF16.h"
#include "Leopard2Backend.h"
#include "leopard2.h"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <new>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

#if defined(__linux__)
#include <sched.h>
#endif

#ifndef __has_feature
#define __has_feature(feature) 0
#endif

#if __has_feature(address_sanitizer) || defined(__SANITIZE_ADDRESS__)
#define LEO2_C7_ASAN_COMPILED 1
#else
#define LEO2_C7_ASAN_COMPILED 0
#endif

#if __has_feature(undefined_behavior_sanitizer)
#define LEO2_C7_UBSAN_COMPILED 1
#else
#define LEO2_C7_UBSAN_COMPILED 0
#endif

#if defined(LEO2_C7_REQUIRE_ASAN_UBSAN)
#if !LEO2_C7_ASAN_COMPILED
#error "LEO2_C7_REQUIRE_ASAN_UBSAN requires compiler ASan instrumentation"
#endif
#if !LEO2_C7_UBSAN_COMPILED
#error "LEO2_C7_REQUIRE_ASAN_UBSAN requires compiler UBSan instrumentation"
#endif
#endif

static bool C7TrackAllocations = false;
static uint64_t C7TrackedAllocations = 0;

#if !defined(LEO2_C7_DISABLE_GLOBAL_NEW_TRACKING)
#if defined(__GNUC__) || defined(__clang__)
#define C7_NOINLINE __attribute__((noinline))
#else
#define C7_NOINLINE
#endif
static C7_NOINLINE void* C7Allocate(size_t bytes)
{
    return malloc(bytes ? bytes : 1);
}
static C7_NOINLINE void C7Deallocate(void* pointer)
{
    free(pointer);
}
void* operator new(size_t bytes)
{
    if (C7TrackAllocations)
        ++C7TrackedAllocations;
    void* pointer = C7Allocate(bytes);
    if (!pointer)
        throw std::bad_alloc();
    return pointer;
}
void* operator new[](size_t bytes) { return ::operator new(bytes); }
void operator delete(void* pointer) noexcept { C7Deallocate(pointer); }
void operator delete[](void* pointer) noexcept { C7Deallocate(pointer); }
void operator delete(void* pointer, size_t) noexcept { C7Deallocate(pointer); }
void operator delete[](void* pointer, size_t) noexcept { C7Deallocate(pointer); }
#define LEO2_C7_ALLOCATION_TRACKING_MODE "global-new"
#else
#define LEO2_C7_ALLOCATION_TRACKING_MODE "disabled-for-sanitizer"
#endif

#ifndef LEO2_C7_SOURCE_SHA256
#define LEO2_C7_SOURCE_SHA256 "unbound"
#endif
#ifndef LEO2_C7_CORE_GIT_SHA
#define LEO2_C7_CORE_GIT_SHA "unbound"
#endif
#ifndef LEO2_C7_LIBRARY_SHA256
#define LEO2_C7_LIBRARY_SHA256 "unbound"
#endif
#ifndef LEO2_C7_SANITIZER_MODE
#define LEO2_C7_SANITIZER_MODE "none"
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

static leo2_backend ParseBackend(const std::string& name)
{
    if (name == "auto")
        return LEO2_BACKEND_AUTO;
    if (name == "scalar")
        return LEO2_BACKEND_SCALAR;
    if (name == "ssse3")
        return LEO2_BACKEND_SSSE3;
    if (name == "avx2")
        return LEO2_BACKEND_AVX2;
    Fail("unknown C7 backend name");
    return LEO2_BACKEND_AUTO;
}

static std::vector<unsigned> ProcessAffinity()
{
    std::vector<unsigned> result;
#if defined(__linux__)
    cpu_set_t set;
    CPU_ZERO(&set);
    if (sched_getaffinity(0, sizeof(set), &set) != 0)
        Fail("sched_getaffinity failed");
    for (unsigned cpu = 0; cpu < CPU_SETSIZE; ++cpu)
        if (CPU_ISSET(cpu, &set))
            result.push_back(cpu);
#endif
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

static bool RangesOverlap(const void* left, const void* right, uint64_t bytes)
{
    const uintptr_t a = reinterpret_cast<uintptr_t>(left);
    const uintptr_t b = reinterpret_cast<uintptr_t>(right);
    if (bytes > std::numeric_limits<uintptr_t>::max() - a ||
        bytes > std::numeric_limits<uintptr_t>::max() - b)
        return true;
    return a < b + static_cast<uintptr_t>(bytes) &&
           b < a + static_cast<uintptr_t>(bytes);
}

template <typename Element> struct FieldOps;

template <> struct FieldOps<uint8_t>
{
    typedef uint8_t Element;
    static unsigned Order() { return leopard::ff8::kOrder; }
    static unsigned Bits() { return 8; }
    static const char* Name() { return "gf8"; }
    static Element Multiply(Element a, Element b)
        { return leopard::ff8::MultiplyElements(a, b); }
    static Element Inverse(Element value)
        { return leopard::ff8::InverseElement(value); }
    static Element Log(Element value)
        { return leopard::ff8::ElementLog(value); }
    static Element FromLog(Element log)
        { return leopard::ff8::MultiplyLogElement(1, log); }
    static void MultiplyBytes(const leopard::backend::Ops& backend_ops,
                              void* output, const void* input, Element log,
                              uint64_t bytes)
        { leopard::ff8::MultiplyBytes(
              backend_ops, output, input, log, bytes); }
    static void MultiplyAddBytes(const leopard::backend::Ops& backend_ops,
                                 void* output, const void* input, Element log,
                                 uint64_t bytes)
        { leopard::ff8::MultiplyAddBytes(
              backend_ops, output, input, log, bytes); }
};

template <> struct FieldOps<uint16_t>
{
    typedef uint16_t Element;
    static unsigned Order() { return leopard::ff16::kOrder; }
    static unsigned Bits() { return 16; }
    static const char* Name() { return "gf16"; }
    static Element Multiply(Element a, Element b)
        { return leopard::ff16::MultiplyElements(a, b); }
    static Element Inverse(Element value)
        { return leopard::ff16::InverseElement(value); }
    static Element Log(Element value)
        { return leopard::ff16::ElementLog(value); }
    static Element FromLog(Element log)
        { return leopard::ff16::MultiplyLogElement(1, log); }
    static void MultiplyBytes(const leopard::backend::Ops& backend_ops,
                              void* output, const void* input, Element log,
                              uint64_t bytes)
        { leopard::ff16::MultiplyBytes(
              backend_ops, output, input, log, bytes); }
    static void MultiplyAddBytes(const leopard::backend::Ops& backend_ops,
                                 void* output, const void* input, Element log,
                                 uint64_t bytes)
        { leopard::ff16::MultiplyAddBytes(
              backend_ops, output, input, log, bytes); }
};

template <typename Element>
class IndependentField
{
public:
    IndependentField(unsigned bits, unsigned polynomial,
                     const Element* basis)
        : bits_(bits), polynomial_(polynomial), order_(1U << bits)
        , to_polynomial_(order_), to_coordinate_(order_)
    {
        std::fill(to_coordinate_.begin(), to_coordinate_.end(),
                  static_cast<Element>(order_ - 1));
        std::vector<uint8_t> seen(order_, 0);
        for (unsigned coordinate = 0; coordinate < order_; ++coordinate)
        {
            unsigned polynomial_value = 0;
            for (unsigned bit = 0; bit < bits_; ++bit)
                if (coordinate & (1U << bit))
                    polynomial_value ^= basis[bit];
            if (polynomial_value >= order_ || seen[polynomial_value])
                Fail("independent coordinate basis is not bijective");
            seen[polynomial_value] = 1;
            to_polynomial_[coordinate] = static_cast<Element>(polynomial_value);
            to_coordinate_[polynomial_value] = static_cast<Element>(coordinate);
        }
    }

    Element Multiply(Element left, Element right) const
    {
        uint32_t a = to_polynomial_[left];
        uint32_t b = to_polynomial_[right];
        uint32_t product = 0;
        while (b)
        {
            if (b & 1U)
                product ^= a;
            b >>= 1;
            a <<= 1;
        }
        for (int bit = static_cast<int>(bits_ * 2 - 2);
             bit >= static_cast<int>(bits_); --bit)
            if (product & (UINT32_C(1) << bit))
                product ^= polynomial_ << (bit - bits_);
        return to_coordinate_[product];
    }

    Element Power(Element value, unsigned exponent) const
    {
        Element result = 1;
        while (exponent)
        {
            if (exponent & 1U)
                result = Multiply(result, value);
            exponent >>= 1;
            if (exponent)
                value = Multiply(value, value);
        }
        return result;
    }

    Element Inverse(Element value) const
    {
        if (!value)
            Fail("independent inverse of zero");
        return Power(value, order_ - 2);
    }

private:
    unsigned bits_;
    uint32_t polynomial_;
    unsigned order_;
    std::vector<Element> to_polynomial_;
    std::vector<Element> to_coordinate_;
};

static IndependentField<uint8_t> IndependentGF8()
{
    static const uint8_t basis[8] = { 1, 214, 152, 146, 86, 200, 88, 230 };
    return IndependentField<uint8_t>(8, 0x11d, basis);
}

static IndependentField<uint16_t> IndependentGF16()
{
    static const uint16_t basis[16] = {
        0x0001, 0xACCA, 0x3C0E, 0x163E,
        0xC582, 0xED2E, 0x914C, 0x4012,
        0x6C98, 0x10D8, 0x6A72, 0xB900,
        0xFDB8, 0xFB34, 0xFF38, 0x991E
    };
    return IndependentField<uint16_t>(16, 0x1002d, basis);
}

template <typename Element>
static std::vector<Element> IndependentRows(
    const IndependentField<Element>& field, unsigned k, unsigned r)
{
    std::vector<Element> weights(k);
    std::vector<Element> rows(static_cast<size_t>(k) * r);
    for (unsigned i = 0; i < k; ++i)
    {
        Element denominator = 1;
        for (unsigned other = 0; other < k; ++other)
            if (other != i)
                denominator = field.Multiply(
                    denominator, static_cast<Element>(i ^ other));
        weights[i] = field.Inverse(denominator);
    }
    for (unsigned parity = 0; parity < r; ++parity)
    {
        const Element point = static_cast<Element>(k + parity);
        Element vanishing = 1;
        for (unsigned i = 0; i < k; ++i)
            vanishing = field.Multiply(
                vanishing, static_cast<Element>(point ^ i));
        for (unsigned i = 0; i < k; ++i)
            rows[static_cast<size_t>(parity) * k + i] = field.Multiply(
                field.Multiply(vanishing,
                    field.Inverse(static_cast<Element>(point ^ i))),
                weights[i]);
    }
    return rows;
}

template <typename Element>
static std::vector<Element> IndependentVandermondeRows(
    const IndependentField<Element>& field, unsigned k, unsigned r)
{
    std::vector<std::vector<Element> > augmented(
        k, std::vector<Element>(static_cast<size_t>(k) * 2, 0));
    for (unsigned row = 0; row < k; ++row)
    {
        Element power = 1;
        for (unsigned column = 0; column < k; ++column)
        {
            augmented[row][column] = power;
            power = field.Multiply(power, static_cast<Element>(row));
        }
        augmented[row][k + row] = 1;
    }
    for (unsigned column = 0; column < k; ++column)
    {
        unsigned pivot = column;
        while (pivot < k && augmented[pivot][column] == 0)
            ++pivot;
        if (pivot == k)
            Fail("independent Vandermonde matrix is singular");
        if (pivot != column)
            augmented[pivot].swap(augmented[column]);
        const Element inverse = field.Inverse(augmented[column][column]);
        for (unsigned value = 0; value < k * 2; ++value)
            augmented[column][value] = field.Multiply(
                augmented[column][value], inverse);
        for (unsigned row = 0; row < k; ++row)
        {
            if (row == column || augmented[row][column] == 0)
                continue;
            const Element factor = augmented[row][column];
            for (unsigned value = 0; value < k * 2; ++value)
                augmented[row][value] ^= field.Multiply(
                    factor, augmented[column][value]);
        }
    }
    std::vector<Element> rows(static_cast<size_t>(k) * r, 0);
    for (unsigned parity = 0; parity < r; ++parity)
    {
        const Element point = static_cast<Element>(k + parity);
        std::vector<Element> powers(k, 1);
        for (unsigned degree = 1; degree < k; ++degree)
            powers[degree] = field.Multiply(powers[degree - 1], point);
        for (unsigned original = 0; original < k; ++original)
            for (unsigned degree = 0; degree < k; ++degree)
                rows[static_cast<size_t>(parity) * k + original] ^=
                    field.Multiply(powers[degree],
                                   augmented[degree][k + original]);
    }
    return rows;
}

template <typename Ops>
struct ExactCodec
{
    typedef typename Ops::Element Element;
    unsigned k;
    unsigned r;
    const leopard::backend::Ops* backend_ops;
    std::vector<Element> logs;

    ExactCodec(const leopard::backend::Ops& selected_backend,
               unsigned original_count, unsigned recovery_count)
        : k(original_count), r(recovery_count)
        , backend_ops(&selected_backend)
        , logs(static_cast<size_t>(k) * r)
    {
        if (k == 0 || r == 0 || static_cast<uint64_t>(k) + r > Ops::Order())
            Fail("exact-low counts exceed the selected field");
        std::vector<Element> weights(k);
        for (unsigned i = 0; i < k; ++i)
        {
            Element denominator = 1;
            for (unsigned other = 0; other < k; ++other)
                if (other != i)
                    denominator = Ops::Multiply(
                        denominator, static_cast<Element>(i ^ other));
            weights[i] = Ops::Inverse(denominator);
        }
        for (unsigned parity = 0; parity < r; ++parity)
        {
            const Element point = static_cast<Element>(k + parity);
            Element vanishing = 1;
            for (unsigned i = 0; i < k; ++i)
                vanishing = Ops::Multiply(
                    vanishing, static_cast<Element>(point ^ i));
            for (unsigned i = 0; i < k; ++i)
            {
                Element coefficient = Ops::Multiply(
                    Ops::Multiply(vanishing,
                        Ops::Inverse(static_cast<Element>(point ^ i))),
                    weights[i]);
                if (!coefficient)
                    Fail("exact-low generator unexpectedly contains zero");
                logs[static_cast<size_t>(parity) * k + i] = Ops::Log(coefficient);
            }
        }
    }

    Element Coefficient(unsigned parity, unsigned original) const
    {
        return Ops::FromLog(logs[static_cast<size_t>(parity) * k + original]);
    }

    // Null parity outputs implement an exact requested subset.  All setup and
    // branch choices are outside the byte inner loops.
    void Encode(uint64_t bytes, const void* const* original,
                void* const* recovery) const
    {
        if (bytes == 0 || bytes > std::numeric_limits<size_t>::max() ||
            (sizeof(Element) == 2 && (bytes & 1U)) ||
            !original || !recovery)
            Fail("invalid exact-low encode arguments or GF16 symbol tail");
        for (unsigned i = 0; i < k; ++i)
            if (!original[i])
                Fail("exact-low source pointer is null");
        for (unsigned parity = 0; parity < r; ++parity)
        {
            if (!recovery[parity])
                continue;
            for (unsigned original_index = 0; original_index < k; ++original_index)
                if (RangesOverlap(recovery[parity], original[original_index], bytes))
                    Fail("exact-low encode output overlaps an input");
            for (unsigned earlier = 0; earlier < parity; ++earlier)
                if (recovery[earlier] &&
                    RangesOverlap(recovery[parity], recovery[earlier], bytes))
                    Fail("exact-low parity outputs overlap");
        }
        for (unsigned parity = 0; parity < r; ++parity)
        {
            void* output = recovery[parity];
            if (!output)
                continue;
            const Element first = logs[static_cast<size_t>(parity) * k];
            if (first == 0)
                memcpy(output, original[0], static_cast<size_t>(bytes));
            else
                Ops::MultiplyBytes(
                    *backend_ops, output, original[0], first, bytes);
            for (unsigned i = 1; i < k; ++i)
                Ops::MultiplyAddBytes(
                    *backend_ops, output, original[i],
                    logs[static_cast<size_t>(parity) * k + i], bytes);
        }
    }
};

template <typename Element>
struct DecodeTerm
{
    bool parity;
    uint32_t index;
    Element multiplier_log;
};

template <typename Ops>
struct ExactDecodePlan
{
    typedef typename Ops::Element Element;
    unsigned k;
    unsigned r;
    const leopard::backend::Ops* backend_ops;
    std::vector<unsigned> missing;
    std::vector<unsigned> selected_parities;
    std::vector<size_t> offsets;
    std::vector<DecodeTerm<Element> > terms;

    ExactDecodePlan(const ExactCodec<Ops>& exact_codec,
                    const std::vector<unsigned>& missing_originals,
                    const std::vector<uint8_t>& parity_present)
        : k(exact_codec.k), r(exact_codec.r)
        , backend_ops(exact_codec.backend_ops), missing(missing_originals)
    {
        if (!std::is_sorted(missing.begin(), missing.end()) ||
            std::adjacent_find(missing.begin(), missing.end()) != missing.end() ||
            (!missing.empty() && missing.back() >= k) ||
            missing.size() > r || parity_present.size() != r)
            Fail("invalid exact-low decode pattern");
        offsets.push_back(0);
        if (missing.empty())
            return;
        for (unsigned parity = 0;
             parity < r && selected_parities.size() < missing.size();
             ++parity)
            if (parity_present[parity])
                selected_parities.push_back(parity);
        if (selected_parities.size() != missing.size())
            Fail("not enough exact-low parity equations");

        const size_t losses = missing.size();
        std::vector<std::vector<Element> > augmented(
            losses, std::vector<Element>(losses * 2, 0));
        for (size_t equation = 0; equation < losses; ++equation)
        {
            for (size_t column = 0; column < losses; ++column)
                augmented[equation][column] = exact_codec.Coefficient(
                    selected_parities[equation], missing[column]);
            augmented[equation][losses + equation] = 1;
        }
        for (size_t column = 0; column < losses; ++column)
        {
            size_t pivot = column;
            while (pivot < losses && augmented[pivot][column] == 0)
                ++pivot;
            if (pivot == losses)
                Fail("exact-low repair minor is singular");
            if (pivot != column)
                augmented[pivot].swap(augmented[column]);
            const Element inverse = Ops::Inverse(augmented[column][column]);
            for (size_t value = 0; value < losses * 2; ++value)
                augmented[column][value] = Ops::Multiply(
                    augmented[column][value], inverse);
            for (size_t row = 0; row < losses; ++row)
            {
                if (row == column || augmented[row][column] == 0)
                    continue;
                const Element factor = augmented[row][column];
                for (size_t value = 0; value < losses * 2; ++value)
                    augmented[row][value] ^= Ops::Multiply(
                        factor, augmented[column][value]);
            }
        }

        for (size_t output = 0; output < losses; ++output)
        {
            for (size_t equation = 0; equation < losses; ++equation)
            {
                const Element coefficient = augmented[output][losses + equation];
                if (coefficient)
                {
                    DecodeTerm<Element> term = {
                        true, selected_parities[equation], Ops::Log(coefficient)
                    };
                    terms.push_back(term);
                }
            }
            for (unsigned original = 0; original < k; ++original)
            {
                if (std::binary_search(missing.begin(), missing.end(), original))
                    continue;
                Element coefficient = 0;
                for (size_t equation = 0; equation < losses; ++equation)
                    coefficient ^= Ops::Multiply(
                        augmented[output][losses + equation],
                        exact_codec.Coefficient(
                            selected_parities[equation], original));
                if (coefficient)
                {
                    DecodeTerm<Element> term = {
                        false, original, Ops::Log(coefficient)
                    };
                    terms.push_back(term);
                }
            }
            if (terms.size() == offsets.back())
                Fail("exact-low repair output contains no terms");
            offsets.push_back(terms.size());
        }
    }

    void Execute(uint64_t bytes, const void* const* original,
                 const void* const* recovery, void* const* restored) const
    {
        if (missing.empty())
            return;
        if (bytes == 0 || bytes > std::numeric_limits<size_t>::max() ||
            (sizeof(Element) == 2 && (bytes & 1U)) ||
            !original || !recovery || !restored)
            Fail("invalid exact-low decode arguments or GF16 symbol tail");
        for (size_t output = 0; output < missing.size(); ++output)
        {
            void* destination = restored[missing[output]];
            if (!destination)
                Fail("exact-low restored pointer is null");
            for (size_t earlier = 0; earlier < output; ++earlier)
                if (RangesOverlap(destination, restored[missing[earlier]], bytes))
                    Fail("exact-low restored outputs overlap");
            for (unsigned source = 0; source < k; ++source)
                if (original[source] && RangesOverlap(destination, original[source], bytes))
                    Fail("exact-low restored output overlaps an original input");
            for (unsigned parity = 0; parity < r; ++parity)
                if (recovery[parity] && RangesOverlap(destination, recovery[parity], bytes))
                    Fail("exact-low restored output overlaps a parity input");
        }
        // Validate every input term for every output before writing the first
        // byte.  This makes pointer/precondition failures atomic across a
        // multi-output repair: a late null survivor or selected parity cannot
        // leave earlier restored shards partially committed.
        for (size_t output = 0; output < missing.size(); ++output)
            for (size_t index = offsets[output]; index < offsets[output + 1]; ++index)
            {
                const DecodeTerm<Element>& term = terms[index];
                const void* source = term.parity
                    ? recovery[term.index] : original[term.index];
                if (!source)
                    Fail(term.parity
                        ? "selected exact-low parity pointer is null"
                        : "surviving exact-low original pointer is null");
            }
        for (size_t output = 0; output < missing.size(); ++output)
        {
            void* destination = restored[missing[output]];
            const size_t begin = offsets[output];
            const size_t end = offsets[output + 1];
            for (size_t index = begin; index < end; ++index)
            {
                const DecodeTerm<Element>& term = terms[index];
                const void* source = term.parity
                    ? recovery[term.index] : original[term.index];
                if (index == begin)
                {
                    if (term.multiplier_log == 0)
                        memcpy(destination, source, static_cast<size_t>(bytes));
                    else
                        Ops::MultiplyBytes(
                            *backend_ops, destination, source,
                            term.multiplier_log, bytes);
                }
                else
                    Ops::MultiplyAddBytes(
                        *backend_ops, destination, source,
                        term.multiplier_log, bytes);
            }
        }
    }
};

static const leopard::backend::Ops* C7TraceTargets[3] = { NULL, NULL, NULL };
static std::atomic<uint64_t> C7TraceCalls[3];

template <unsigned Slot>
static void TraceFF8Multiply(
    void* output, const void* input, uint16_t log, uint64_t bytes)
{
    C7TraceCalls[Slot].fetch_add(1, std::memory_order_relaxed);
    C7TraceTargets[Slot]->ff8_multiply(output, input, log, bytes);
}

template <unsigned Slot>
static void TraceFF8MultiplyAdd(
    void* output, const void* input, uint16_t log, uint64_t bytes)
{
    C7TraceCalls[Slot].fetch_add(1, std::memory_order_relaxed);
    C7TraceTargets[Slot]->ff8_multiply_add(output, input, log, bytes);
}

template <unsigned Slot>
static void TraceFF16Multiply(
    void* output, const void* input, uint16_t log, uint64_t bytes)
{
    C7TraceCalls[Slot].fetch_add(1, std::memory_order_relaxed);
    C7TraceTargets[Slot]->ff16_multiply(output, input, log, bytes);
}

template <unsigned Slot>
static void TraceFF16MultiplyAdd(
    void* output, const void* input, uint16_t log, uint64_t bytes)
{
    C7TraceCalls[Slot].fetch_add(1, std::memory_order_relaxed);
    C7TraceTargets[Slot]->ff16_multiply_add(output, input, log, bytes);
}

template <unsigned Slot>
static void ConfigureTracingOps(
    const leopard::backend::Ops& target,
    leopard::backend::Ops& tracing)
{
    C7TraceTargets[Slot] = &target;
    C7TraceCalls[Slot].store(0, std::memory_order_relaxed);
    tracing = target;
    tracing.ff8_multiply = &TraceFF8Multiply<Slot>;
    tracing.ff8_multiply_add = &TraceFF8MultiplyAdd<Slot>;
    tracing.ff16_multiply = &TraceFF16Multiply<Slot>;
    tracing.ff16_multiply_add = &TraceFF16MultiplyAdd<Slot>;
}

template <typename Element>
static Element LoadSymbol(const uint8_t* data, size_t bytes, size_t symbol)
{
    if (sizeof(Element) == 1)
        return static_cast<Element>(data[symbol]);
    const size_t complete = bytes & ~static_cast<size_t>(63);
    const size_t complete_symbols = complete / 2;
    if (symbol < complete_symbols)
    {
        const size_t tile = (symbol / 32) * 64;
        const size_t lane = symbol % 32;
        return static_cast<Element>(data[tile + lane] |
            (static_cast<unsigned>(data[tile + 32 + lane]) << 8));
    }
    const size_t residual_symbols = (bytes - complete) / 2;
    const size_t lane = symbol - complete_symbols;
    if (lane >= residual_symbols)
        Fail("GF16 symbol index exceeds compact tail");
    return static_cast<Element>(data[complete + lane] |
        (static_cast<unsigned>(data[complete + residual_symbols + lane]) << 8));
}

template <typename Element>
static void StoreSymbol(uint8_t* data, size_t bytes, size_t symbol, Element value)
{
    if (sizeof(Element) == 1)
    {
        data[symbol] = static_cast<uint8_t>(value);
        return;
    }
    const size_t complete = bytes & ~static_cast<size_t>(63);
    const size_t complete_symbols = complete / 2;
    if (symbol < complete_symbols)
    {
        const size_t tile = (symbol / 32) * 64;
        const size_t lane = symbol % 32;
        data[tile + lane] = static_cast<uint8_t>(value);
        data[tile + 32 + lane] = static_cast<uint8_t>(value >> 8);
        return;
    }
    const size_t residual_symbols = (bytes - complete) / 2;
    const size_t lane = symbol - complete_symbols;
    if (lane >= residual_symbols)
        Fail("GF16 symbol index exceeds compact tail");
    data[complete + lane] = static_cast<uint8_t>(value);
    data[complete + residual_symbols + lane] = static_cast<uint8_t>(value >> 8);
}

template <typename Ops>
static uint64_t RunConcurrentBackendCase(
    const leopard::backend::Ops& backend_ops, size_t bytes,
    uint64_t& executions)
{
    const unsigned k = 8;
    const unsigned r = 8;
    const unsigned iterations = 32;
    ExactCodec<Ops> codec(backend_ops, k, r);
    std::vector<std::vector<uint8_t> > source(
        k, std::vector<uint8_t>(bytes));
    std::vector<std::vector<uint8_t> > parity(
        r, std::vector<uint8_t>(bytes));
    std::vector<std::vector<uint8_t> > restored(
        k, std::vector<uint8_t>(bytes, 0xcc));
    std::vector<const void*> source_pointer(k);
    std::vector<void*> parity_pointer(r);
    for (unsigned original = 0; original < k; ++original)
    {
        for (size_t byte = 0; byte < bytes; ++byte)
            source[original][byte] = static_cast<uint8_t>(
                original * 73U + static_cast<unsigned>(byte) * 29U + 11U);
        source_pointer[original] = &source[original][0];
    }
    for (unsigned recovery = 0; recovery < r; ++recovery)
        parity_pointer[recovery] = &parity[recovery][0];

    const unsigned missing_array[] = { 1, 6 };
    const std::vector<unsigned> missing(
        missing_array, missing_array + sizeof(missing_array) / sizeof(unsigned));
    std::vector<uint8_t> parity_present(r, 1);
    parity_present[0] = 0;
    const ExactDecodePlan<Ops> plan(codec, missing, parity_present);
    std::vector<const void*> received_source = source_pointer;
    received_source[missing[0]] = NULL;
    received_source[missing[1]] = NULL;
    std::vector<const void*> received_parity(r);
    for (unsigned recovery = 0; recovery < r; ++recovery)
        received_parity[recovery] = parity_present[recovery]
            ? static_cast<const void*>(&parity[recovery][0]) : NULL;
    std::vector<void*> restored_pointer(k, NULL);
    restored_pointer[missing[0]] = &restored[missing[0]][0];
    restored_pointer[missing[1]] = &restored[missing[1]][0];

    for (unsigned iteration = 0; iteration < iterations; ++iteration)
    {
        codec.Encode(bytes, &source_pointer[0], &parity_pointer[0]);
        plan.Execute(bytes, &received_source[0], &received_parity[0],
                     &restored_pointer[0]);
        for (size_t i = 0; i < missing.size(); ++i)
            if (memcmp(&restored[missing[i]][0], &source[missing[i]][0],
                       bytes) != 0)
                Fail("concurrent exact-low plan restored the wrong original");
        executions += 2;
    }

    uint64_t digest = UINT64_C(1469598103934665603);
    for (unsigned recovery = 0; recovery < r; ++recovery)
        digest = Fnv(digest, &parity[recovery][0], bytes);
    for (size_t i = 0; i < missing.size(); ++i)
        digest = Fnv(digest, &restored[missing[i]][0], bytes);
    return digest;
}

struct Correctness
{
    uint64_t gf8_cases;
    uint64_t gf16_cases;
    uint64_t coefficients;
    uint64_t gf16_vandermonde_coefficients;
    uint64_t encode_executions;
    uint64_t encode_symbol_comparisons;
    uint64_t subset_encode_executions;
    uint64_t decode_plans;
    uint64_t decode_executions;
    uint64_t decode_symbol_comparisons;
    uint64_t maximum_loss_plans;
    uint64_t unavailable_parity_plans;
    uint64_t no_loss_null_calls;
    uint64_t parity_rebuilds;
    uint64_t odd_gf16_rejections;
    uint64_t overlap_rejections;
    uint64_t parity_output_overlap_rejections;
    uint64_t restored_output_overlap_rejections;
    uint64_t restored_input_overlap_rejections;
    uint64_t selected_parity_null_rejections;
    uint64_t survivor_null_rejections;
    uint64_t atomic_rejection_bytes_checked;
    uint64_t read_only_input_alias_calls;
    uint64_t read_only_input_alias_symbol_comparisons;
    uint64_t decode_read_only_input_alias_calls;
    uint64_t decode_read_only_input_alias_symbol_comparisons;
    uint64_t detached_plan_executions;
    uint64_t detached_plan_symbol_comparisons;
    uint64_t concurrent_backend_contexts;
    uint64_t concurrent_backend_executions;
    uint64_t concurrent_backend_trace_calls;
    uint64_t concurrent_wire_digest_comparisons;
    uint64_t exhaustive_small_plans;
    uint64_t exhaustive_small_executions;
    uint64_t exhaustive_small_symbol_comparisons;
    uint64_t malformed_plan_rejections;
    uint64_t hot_path_allocations;
    uint64_t digest;
};

static void ValidateConcurrentBackendContexts(Correctness& result)
{
    const leo2_backend kinds[3] = {
        LEO2_BACKEND_SCALAR, LEO2_BACKEND_SSSE3, LEO2_BACKEND_AVX2
    };
    leo2_context* contexts[3] = { NULL, NULL, NULL };
    leopard::backend::Ops tracing[3] = {};
    bool available[3] = { false, false, false };
    for (unsigned slot = 0; slot < 3; ++slot)
    {
        leo2_context_options options = {};
        options.struct_size = sizeof(options);
        options.backend = static_cast<uint32_t>(kinds[slot]);
        options.thread_count = 1;
        const leo2_result created = leo2_context_create(&options, &contexts[slot]);
        if (created == LEO2_UNSUPPORTED)
            continue;
        Require(created, "C7 concurrent backend context create");
        if (leo2_context_backend(contexts[slot]) != kinds[slot])
            Fail("concurrent C7 context selected the wrong backend");
        leopard::backend::QualificationStatus qualification =
            leopard::backend::QualificationAvailable;
        const leopard::backend::Ops* target =
            leopard::backend::GetQualifiedOps(kinds[slot], &qualification);
        if (!target || target->kind != kinds[slot])
            Fail("concurrent C7 backend table did not qualify exactly");
        if (slot == 0)
            ConfigureTracingOps<0>(*target, tracing[slot]);
        else if (slot == 1)
            ConfigureTracingOps<1>(*target, tracing[slot]);
        else
            ConfigureTracingOps<2>(*target, tracing[slot]);
        available[slot] = true;
        ++result.concurrent_backend_contexts;
    }
    if (!available[0])
        Fail("scalar concurrent C7 backend is unavailable");

    uint64_t gf8_digest[3] = {};
    uint64_t gf16_digest[3] = {};
    uint64_t executions[3] = {};
    std::exception_ptr errors[3];
    std::vector<std::thread> threads;
    for (unsigned slot = 0; slot < 3; ++slot)
    {
        if (!available[slot])
            continue;
        threads.push_back(std::thread([&, slot]() {
            try
            {
                gf8_digest[slot] = RunConcurrentBackendCase<FieldOps<uint8_t> >(
                    tracing[slot], 129, executions[slot]);
                gf16_digest[slot] = RunConcurrentBackendCase<FieldOps<uint16_t> >(
                    tracing[slot], 130, executions[slot]);
            }
            catch (...)
            {
                errors[slot] = std::current_exception();
            }
        }));
    }
    for (size_t thread = 0; thread < threads.size(); ++thread)
        threads[thread].join();

    unsigned reference = 0;
    while (reference < 3 && !available[reference])
        ++reference;
    for (unsigned slot = 0; slot < 3; ++slot)
    {
        if (!available[slot])
            continue;
        if (errors[slot])
            std::rethrow_exception(errors[slot]);
        if (gf8_digest[slot] != gf8_digest[reference] ||
            gf16_digest[slot] != gf16_digest[reference])
            Fail("concurrent exact-low backends changed the wire result");
        if (slot != reference)
            result.concurrent_wire_digest_comparisons += 2;
        const uint64_t calls = C7TraceCalls[slot].load(
            std::memory_order_relaxed);
        if (calls == 0)
            Fail("exact-low execution bypassed its selected backend table");
        result.concurrent_backend_trace_calls += calls;
        result.concurrent_backend_executions += executions[slot];
        leo2_context_destroy(contexts[slot]);
        contexts[slot] = NULL;
    }
}

static unsigned CountBits(unsigned value)
{
    unsigned count = 0;
    while (value)
    {
        count += value & 1U;
        value >>= 1;
    }
    return count;
}

static void ValidateExhaustiveSmallPlans(
    const leopard::backend::Ops& backend_ops, Correctness& result)
{
    typedef FieldOps<uint8_t> Ops;
    const unsigned k = 4;
    const unsigned r = 4;
    const uint64_t bytes = 1;
    ExactCodec<Ops> codec(backend_ops, k, r);
    std::vector<std::vector<uint8_t> > source(k, std::vector<uint8_t>(1));
    std::vector<std::vector<uint8_t> > parity(r, std::vector<uint8_t>(1));
    std::vector<const void*> source_pointer(k);
    std::vector<void*> parity_pointer(r);
    for (unsigned original = 0; original < k; ++original)
        source_pointer[original] = &source[original][0];
    for (unsigned recovery = 0; recovery < r; ++recovery)
        parity_pointer[recovery] = &parity[recovery][0];

    for (unsigned missing_mask = 0; missing_mask < (1U << k); ++missing_mask)
    {
        const unsigned losses = CountBits(missing_mask);
        std::vector<unsigned> missing;
        for (unsigned original = 0; original < k; ++original)
            if (missing_mask & (1U << original))
                missing.push_back(original);
        for (unsigned parity_mask = 0; parity_mask < (1U << r); ++parity_mask)
        {
            if (CountBits(parity_mask) < losses)
                continue;
            std::vector<uint8_t> parity_present(r, 0);
            for (unsigned recovery = 0; recovery < r; ++recovery)
                parity_present[recovery] = static_cast<uint8_t>(
                    (parity_mask >> recovery) & 1U);
            const ExactDecodePlan<Ops> plan(codec, missing, parity_present);
            ++result.exhaustive_small_plans;
            if (missing.empty())
            {
                plan.Execute(bytes, NULL, NULL, NULL);
                ++result.exhaustive_small_executions;
                continue;
            }
            for (unsigned basis = 0; basis < k * 8; ++basis)
            {
                for (unsigned original = 0; original < k; ++original)
                    source[original][0] = 0;
                source[basis / 8][0] = static_cast<uint8_t>(1U << (basis & 7U));
                codec.Encode(bytes, &source_pointer[0], &parity_pointer[0]);
                std::vector<const void*> received_source = source_pointer;
                std::vector<const void*> received_parity(r, NULL);
                std::vector<std::vector<uint8_t> > restored(
                    k, std::vector<uint8_t>(1, 0xa7));
                std::vector<void*> restored_pointer(k, NULL);
                for (size_t i = 0; i < missing.size(); ++i)
                {
                    received_source[missing[i]] = NULL;
                    restored_pointer[missing[i]] = &restored[missing[i]][0];
                }
                for (unsigned recovery = 0; recovery < r; ++recovery)
                    if (parity_present[recovery])
                        received_parity[recovery] = &parity[recovery][0];
                plan.Execute(bytes, &received_source[0], &received_parity[0],
                             &restored_pointer[0]);
                ++result.exhaustive_small_executions;
                for (size_t i = 0; i < missing.size(); ++i)
                {
                    if (restored[missing[i]][0] != source[missing[i]][0])
                        Fail("exhaustive exact-low plan restored wrong basis data");
                    ++result.exhaustive_small_symbol_comparisons;
                }
            }
        }
    }
}

static void RequireFilled(
    const std::vector<std::vector<uint8_t> >& storage,
    uint8_t expected, uint64_t& checked)
{
    for (size_t shard = 0; shard < storage.size(); ++shard)
        for (size_t byte = 0; byte < storage[shard].size(); ++byte)
        {
            if (storage[shard][byte] != expected)
                Fail("rejected execution partially changed an output");
            ++checked;
        }
}

template <typename Callable>
static void ExpectRejected(Callable callable, uint64_t& counter)
{
    try
    {
        callable();
    }
    catch (const std::runtime_error&)
    {
        ++counter;
        return;
    }
    Fail("invalid exact-low execution was not rejected");
}

static std::vector<unsigned> MissingSet(unsigned k, unsigned losses)
{
    std::vector<unsigned> missing;
    for (unsigned i = 0; i < losses; ++i)
        missing.push_back(static_cast<unsigned>(
            (static_cast<uint64_t>(i) * k) / losses));
    if (!std::is_sorted(missing.begin(), missing.end()) ||
        std::adjacent_find(missing.begin(), missing.end()) != missing.end())
        Fail("deterministic missing set is not unique");
    return missing;
}

static uint64_t NextPatternValue(uint64_t& state)
{
    state ^= state >> 12;
    state ^= state << 25;
    state ^= state >> 27;
    return state * UINT64_C(2685821657736338717);
}

static std::vector<unsigned> VariedMissingSet(
    unsigned k, unsigned losses, uint64_t seed)
{
    if (losses == 0 || losses > k)
        Fail("invalid varied exact-low loss count");
    std::vector<unsigned> coordinates(k);
    std::iota(coordinates.begin(), coordinates.end(), 0U);
    uint64_t state = seed ? seed : UINT64_C(0x9e3779b97f4a7c15);
    for (unsigned i = 0; i < losses; ++i)
    {
        const unsigned selected = i + static_cast<unsigned>(
            NextPatternValue(state) % (k - i));
        std::swap(coordinates[i], coordinates[selected]);
    }
    coordinates.resize(losses);
    std::sort(coordinates.begin(), coordinates.end());
    return coordinates;
}

static bool MarkVariedUnavailableParities(
    std::vector<uint8_t>& parity_present, unsigned required,
    uint64_t seed)
{
    const unsigned r = static_cast<unsigned>(parity_present.size());
    if (required >= r)
        return false;
    const unsigned cap = std::min(7U, r - required);
    uint64_t state = seed ? seed : UINT64_C(0xd1b54a32d192ed03);
    const uint64_t selector = NextPatternValue(state);
    const unsigned unavailable = (selector & 1U)
        ? r - required
        : 1U + static_cast<unsigned>(selector % cap);
    std::vector<unsigned> coordinates(r);
    std::iota(coordinates.begin(), coordinates.end(), 0U);
    for (unsigned i = 0; i < unavailable; ++i)
    {
        const unsigned selected = i + static_cast<unsigned>(
            NextPatternValue(state) % (r - i));
        std::swap(coordinates[i], coordinates[selected]);
        parity_present[coordinates[i]] = 0;
    }
    return true;
}

static std::vector<unsigned> LossCounts(unsigned k, unsigned r)
{
    const unsigned maximum = std::min(k, r);
    const unsigned candidates[] = {
        1, std::min(2U, maximum), std::min(4U, maximum), maximum
    };
    std::vector<unsigned> result(candidates, candidates + 4);
    std::sort(result.begin(), result.end());
    result.erase(std::unique(result.begin(), result.end()), result.end());
    return result;
}

template <typename Ops>
static void ValidateCase(
    const leopard::backend::Ops& backend_ops,
    unsigned k, unsigned r, unsigned case_index,
    const std::vector<size_t>& byte_counts,
    const IndependentField<typename Ops::Element>& independent,
    Correctness& result)
{
    typedef typename Ops::Element Element;
    ExactCodec<Ops> codec(backend_ops, k, r);
    const std::vector<Element> rows = IndependentRows(independent, k, r);
    if (sizeof(Element) == 2 && k == 3 && r == 500)
    {
        const std::vector<Element> vandermonde =
            IndependentVandermondeRows(independent, k, r);
        if (vandermonde != rows)
            Fail("declared GF16 Vandermonde oracle disagrees with barycentric rows");
        result.gf16_vandermonde_coefficients += vandermonde.size();
    }
    for (unsigned parity = 0; parity < r; ++parity)
        for (unsigned original = 0; original < k; ++original)
        {
            const Element actual = codec.Coefficient(parity, original);
            const Element expected = rows[static_cast<size_t>(parity) * k + original];
            if (!actual || actual != expected)
                Fail("production-field and independent exact rows disagree");
            ++result.coefficients;
        }

    std::vector<uint8_t> all_parity_present(r, 1);
    if (k > 1)
    {
        const unsigned unsorted_array[] = { 1, 0 };
        const std::vector<unsigned> unsorted(
            unsorted_array, unsorted_array + 2);
        ExpectRejected([&]() {
            ExactDecodePlan<Ops> invalid(codec, unsorted, all_parity_present);
        }, result.malformed_plan_rejections);
        const std::vector<unsigned> duplicate(2, 0);
        ExpectRejected([&]() {
            ExactDecodePlan<Ops> invalid(codec, duplicate, all_parity_present);
        }, result.malformed_plan_rejections);
    }
    const std::vector<unsigned> out_of_range(1, k);
    ExpectRejected([&]() {
        ExactDecodePlan<Ops> invalid(codec, out_of_range, all_parity_present);
    }, result.malformed_plan_rejections);
    std::vector<uint8_t> wrong_parity_size(r - 1, 1);
    ExpectRejected([&]() {
        ExactDecodePlan<Ops> invalid(
            codec, std::vector<unsigned>(), wrong_parity_size);
    }, result.malformed_plan_rejections);
    std::vector<uint8_t> no_parity(r, 0);
    const std::vector<unsigned> one_missing(1, 0);
    ExpectRejected([&]() {
        ExactDecodePlan<Ops> invalid(codec, one_missing, no_parity);
    }, result.malformed_plan_rejections);
    if (k > r)
    {
        std::vector<unsigned> too_many(r + 1);
        std::iota(too_many.begin(), too_many.end(), 0U);
        ExpectRejected([&]() {
            ExactDecodePlan<Ops> invalid(codec, too_many, all_parity_present);
        }, result.malformed_plan_rejections);
    }
    const ExactDecodePlan<Ops> no_loss(
        codec, std::vector<unsigned>(), all_parity_present);
    const std::vector<unsigned> detached_missing(
        1, static_cast<unsigned>((case_index * 37U + k / 2U) % k));
    std::unique_ptr<ExactDecodePlan<Ops> > detached_plan;
    {
        ExactCodec<Ops> temporary_codec(backend_ops, k, r);
        detached_plan.reset(new ExactDecodePlan<Ops>(
            temporary_codec, detached_missing, all_parity_present));
    }
    C7TrackedAllocations = 0;
    C7TrackAllocations = true;
    no_loss.Execute(1, NULL, NULL, NULL);
    C7TrackAllocations = false;
    result.hot_path_allocations += C7TrackedAllocations;
    ++result.no_loss_null_calls;

    for (size_t size_index = 0; size_index < byte_counts.size(); ++size_index)
    {
        const size_t bytes = byte_counts[size_index];
        if (sizeof(Element) == 2 && (bytes & 1U))
            Fail("GF16 test length must contain whole symbols");
        const size_t symbols = bytes / sizeof(Element);
        std::vector<std::vector<uint8_t> > source(
            k, std::vector<uint8_t>(bytes + 2, 0xa5));
        std::vector<std::vector<uint8_t> > parity(
            r, std::vector<uint8_t>(bytes + 2, 0xa5));
        std::vector<const void*> source_ptr(k);
        std::vector<void*> parity_ptr(r);
        for (unsigned original = 0; original < k; ++original)
        {
            uint8_t* data = &source[original][1];
            for (size_t symbol = 0; symbol < symbols; ++symbol)
            {
                const uint32_t value = original * 977U +
                    static_cast<unsigned>(symbol) * 131U + case_index * 29U +
                    static_cast<unsigned>(size_index) * 7U + 1U;
                StoreSymbol<Element>(data, bytes, symbol,
                    static_cast<Element>(value % Ops::Order()));
            }
            source_ptr[original] = data;
        }
        for (unsigned recovery = 0; recovery < r; ++recovery)
            parity_ptr[recovery] = &parity[recovery][1];

        if (size_index == 0)
        {
            std::vector<void*> overlapping(r, NULL);
            overlapping[0] = const_cast<void*>(source_ptr[0]);
            ExpectRejected([&]() {
                codec.Encode(bytes, &source_ptr[0], &overlapping[0]);
            }, result.overlap_rejections);
            RequireFilled(parity, 0xa5, result.atomic_rejection_bytes_checked);
            if (r > 1)
            {
                overlapping.assign(r, NULL);
                overlapping[0] = &parity[0][1];
                overlapping[1] = &parity[0][1];
                ExpectRejected([&]() {
                    codec.Encode(bytes, &source_ptr[0], &overlapping[0]);
                }, result.parity_output_overlap_rejections);
                ++result.overlap_rejections;
                RequireFilled(parity, 0xa5, result.atomic_rejection_bytes_checked);
            }
            if (sizeof(Element) == 2)
            {
                ExpectRejected([&]() {
                    codec.Encode(3, &source_ptr[0], &parity_ptr[0]);
                }, result.odd_gf16_rejections);
            }
            if (k > 1)
            {
                std::vector<const void*> aliased_source = source_ptr;
                aliased_source[1] = aliased_source[0];
                std::vector<std::vector<uint8_t> > alias_output(
                    r, std::vector<uint8_t>(bytes + 2, 0x79));
                std::vector<void*> alias_output_ptr(r);
                for (unsigned recovery = 0; recovery < r; ++recovery)
                    alias_output_ptr[recovery] = &alias_output[recovery][1];
                C7TrackedAllocations = 0;
                C7TrackAllocations = true;
                codec.Encode(bytes, &aliased_source[0], &alias_output_ptr[0]);
                C7TrackAllocations = false;
                result.hot_path_allocations += C7TrackedAllocations;
                ++result.read_only_input_alias_calls;
                for (unsigned recovery = 0; recovery < r; ++recovery)
                {
                    if (alias_output[recovery][0] != 0x79 ||
                        alias_output[recovery][bytes + 1] != 0x79)
                        Fail("aliased-input encode changed an output guard");
                    for (size_t symbol = 0; symbol < symbols; ++symbol)
                    {
                        Element expected = 0;
                        for (unsigned original = 0; original < k; ++original)
                            expected ^= independent.Multiply(
                                rows[static_cast<size_t>(recovery) * k + original],
                                LoadSymbol<Element>(
                                    static_cast<const uint8_t*>(aliased_source[original]),
                                    bytes, symbol));
                        if (LoadSymbol<Element>(&alias_output[recovery][1],
                                               bytes, symbol) != expected)
                            Fail("read-only input alias encode result is incorrect");
                        ++result.read_only_input_alias_symbol_comparisons;
                    }
                }
            }
        }

        C7TrackedAllocations = 0;
        C7TrackAllocations = true;
        codec.Encode(bytes, &source_ptr[0], &parity_ptr[0]);
        C7TrackAllocations = false;
        result.hot_path_allocations += C7TrackedAllocations;
        ++result.encode_executions;

        for (unsigned recovery = 0; recovery < r; ++recovery)
        {
            if (parity[recovery][0] != 0xa5 || parity[recovery][bytes + 1] != 0xa5)
                Fail("exact-low encoder changed an output guard");
            for (size_t symbol = 0; symbol < symbols; ++symbol)
            {
                Element expected = 0;
                for (unsigned original = 0; original < k; ++original)
                    expected ^= independent.Multiply(
                        rows[static_cast<size_t>(recovery) * k + original],
                        LoadSymbol<Element>(&source[original][1], bytes, symbol));
                const Element actual = LoadSymbol<Element>(
                    &parity[recovery][1], bytes, symbol);
                if (actual != expected)
                    Fail("exact-low encode differs from independent symbol algebra");
                ++result.encode_symbol_comparisons;
            }
            result.digest = Fnv(result.digest, &parity[recovery][1], bytes);
        }

        // Plans are immutable, self-contained execution objects: destroying
        // the codec used for setup must not invalidate its field/backend table,
        // repair coefficients, or coordinate metadata.
        if (size_index + 1 == byte_counts.size())
        {
            std::vector<const void*> detached_source = source_ptr;
            detached_source[detached_missing[0]] = NULL;
            std::vector<const void*> detached_parity(r);
            for (unsigned recovery = 0; recovery < r; ++recovery)
                detached_parity[recovery] = &parity[recovery][1];
            std::vector<uint8_t> detached_storage(bytes + 2, 0x94);
            std::vector<void*> detached_restored(k, NULL);
            detached_restored[detached_missing[0]] = &detached_storage[1];
            C7TrackedAllocations = 0;
            C7TrackAllocations = true;
            detached_plan->Execute(
                bytes, &detached_source[0], &detached_parity[0],
                &detached_restored[0]);
            C7TrackAllocations = false;
            result.hot_path_allocations += C7TrackedAllocations;
            if (detached_storage.front() != 0x94 ||
                detached_storage.back() != 0x94 ||
                memcmp(&detached_storage[1],
                       &source[detached_missing[0]][1], bytes) != 0)
                Fail("decode plan retained destroyed codec state");
            ++result.detached_plan_executions;
            result.detached_plan_symbol_comparisons += symbols;
        }

        // A nonzero constant-polynomial codeword gives every received shard
        // identical bytes, so all surviving-original and parity inputs can
        // validly share one read-only buffer.  Exercise that strongest
        // input/input alias case while restoring one original into separate
        // guarded storage, and snapshot the complete shared input so an
        // accidental write cannot hide behind unchanged guard bytes.
        {
            const std::vector<unsigned> missing(1, 0);
            const ExactDecodePlan<Ops> alias_plan(
                codec, missing, all_parity_present);
            std::vector<uint8_t> shared_input(bytes + 2, 0);
            shared_input.front() = 0x71;
            shared_input.back() = 0x71;
            for (size_t symbol = 0; symbol < symbols; ++symbol)
            {
                const uint32_t value = 1U +
                    (case_index * 193U + static_cast<unsigned>(size_index) * 17U +
                     static_cast<unsigned>(symbol) * 29U) % (Ops::Order() - 1U);
                StoreSymbol<Element>(
                    &shared_input[1], bytes, symbol,
                    static_cast<Element>(value));
            }
            const std::vector<uint8_t> shared_snapshot = shared_input;
            std::vector<const void*> aliased_original(
                k, static_cast<const void*>(&shared_input[1]));
            std::vector<const void*> aliased_parity(
                r, static_cast<const void*>(&shared_input[1]));
            aliased_original[0] = NULL;
            std::vector<uint8_t> restored_alias(bytes + 2, 0x4d);
            std::vector<void*> restored_alias_ptr(k, NULL);
            restored_alias_ptr[0] = &restored_alias[1];

            C7TrackedAllocations = 0;
            C7TrackAllocations = true;
            alias_plan.Execute(bytes, &aliased_original[0],
                               &aliased_parity[0], &restored_alias_ptr[0]);
            C7TrackAllocations = false;
            result.hot_path_allocations += C7TrackedAllocations;
            ++result.decode_read_only_input_alias_calls;
            if (shared_input != shared_snapshot ||
                restored_alias.front() != 0x4d || restored_alias.back() != 0x4d)
                Fail("decode input alias execution changed shared input or a guard");
            for (size_t symbol = 0; symbol < symbols; ++symbol)
            {
                const Element expected = LoadSymbol<Element>(
                    &shared_snapshot[1], bytes, symbol);
                const Element actual = LoadSymbol<Element>(
                    &restored_alias[1], bytes, symbol);
                if (actual != expected)
                    Fail("decode input alias execution restored wrong constant data");
            }
            result.decode_read_only_input_alias_symbol_comparisons += symbols;
        }

        // Requested parity subset: null outputs remain untouched and selected
        // outputs reproduce the corresponding full-encode bytes exactly.
        std::vector<std::vector<uint8_t> > subset(
            r, std::vector<uint8_t>(bytes + 2, 0x5a));
        std::vector<void*> subset_ptr(r, NULL);
        for (unsigned recovery = 0; recovery < r; ++recovery)
            if ((recovery + case_index) % 3 == 0 || recovery + 1 == r)
                subset_ptr[recovery] = &subset[recovery][1];
        C7TrackedAllocations = 0;
        C7TrackAllocations = true;
        codec.Encode(bytes, &source_ptr[0], &subset_ptr[0]);
        C7TrackAllocations = false;
        result.hot_path_allocations += C7TrackedAllocations;
        ++result.subset_encode_executions;
        for (unsigned recovery = 0; recovery < r; ++recovery)
        {
            if (subset_ptr[recovery])
            {
                if (memcmp(&subset[recovery][1], &parity[recovery][1], bytes) != 0)
                    Fail("requested parity subset differs from full encode");
            }
            else if (std::find(subset[recovery].begin(), subset[recovery].end(),
                               static_cast<uint8_t>(0x5a)) == subset[recovery].end() ||
                     !std::all_of(subset[recovery].begin(), subset[recovery].end(),
                                  [](uint8_t value) { return value == 0x5a; }))
                Fail("unrequested parity output changed");
        }

        const std::vector<unsigned> loss_counts = LossCounts(k, r);
        for (size_t loss_index = 0; loss_index < loss_counts.size(); ++loss_index)
        {
            const unsigned losses = loss_counts[loss_index];
            const uint64_t pattern_seed =
                (static_cast<uint64_t>(case_index + 1U) << 40) ^
                (static_cast<uint64_t>(size_index + 1U) << 20) ^
                static_cast<uint64_t>(loss_index + 1U);
            const std::vector<unsigned> missing = VariedMissingSet(
                k, losses, pattern_seed);
            std::vector<uint8_t> parity_present(r, 1);
            bool has_unavailable_parity = false;
            if (r > losses && ((case_index + size_index + loss_index) & 1U))
                has_unavailable_parity = MarkVariedUnavailableParities(
                    parity_present, losses,
                    pattern_seed ^ UINT64_C(0xa0761d6478bd642f));
            const ExactDecodePlan<Ops> plan(codec, missing, parity_present);
            ++result.decode_plans;
            result.maximum_loss_plans += losses == std::min(k, r);
            result.unavailable_parity_plans += has_unavailable_parity;

            std::vector<const void*> received_source = source_ptr;
            for (size_t i = 0; i < missing.size(); ++i)
                received_source[missing[i]] = NULL;
            std::vector<const void*> received_parity(r);
            for (unsigned recovery = 0; recovery < r; ++recovery)
                received_parity[recovery] = parity_present[recovery]
                    ? static_cast<const void*>(&parity[recovery][1]) : NULL;
            std::vector<std::vector<uint8_t> > restored_storage(
                k, std::vector<uint8_t>(bytes + 2, 0x3c));
            std::vector<void*> restored(k);
            for (unsigned original = 0; original < k; ++original)
                restored[original] = &restored_storage[original][1];

            if (size_index == 0 && loss_index + 1 == loss_counts.size())
            {
                void* saved = restored[missing[0]];
                restored[missing[0]] = const_cast<void*>(
                    received_parity[plan.selected_parities[0]]);
                ExpectRejected([&]() {
                    plan.Execute(bytes, &received_source[0],
                                 &received_parity[0], &restored[0]);
                }, result.restored_input_overlap_rejections);
                ++result.overlap_rejections;
                restored[missing[0]] = saved;
                RequireFilled(restored_storage, 0x3c,
                              result.atomic_rejection_bytes_checked);

                unsigned survivor = codec.k;
                for (unsigned original = 0; original < codec.k; ++original)
                    if (!std::binary_search(missing.begin(), missing.end(), original))
                    {
                        survivor = original;
                        break;
                    }
                if (survivor < codec.k)
                {
                    restored[missing[0]] = const_cast<void*>(
                        received_source[survivor]);
                    ExpectRejected([&]() {
                        plan.Execute(bytes, &received_source[0],
                                     &received_parity[0], &restored[0]);
                    }, result.restored_input_overlap_rejections);
                    ++result.overlap_rejections;
                    restored[missing[0]] = saved;
                    RequireFilled(restored_storage, 0x3c,
                                  result.atomic_rejection_bytes_checked);
                }

                if (missing.size() > 1)
                {
                    void* second = restored[missing[1]];
                    restored[missing[1]] = restored[missing[0]];
                    ExpectRejected([&]() {
                        plan.Execute(bytes, &received_source[0],
                                     &received_parity[0], &restored[0]);
                    }, result.restored_output_overlap_rejections);
                    ++result.overlap_rejections;
                    restored[missing[1]] = second;
                    RequireFilled(restored_storage, 0x3c,
                                  result.atomic_rejection_bytes_checked);
                }

                const unsigned selected = plan.selected_parities.back();
                const void* selected_pointer = received_parity[selected];
                received_parity[selected] = NULL;
                ExpectRejected([&]() {
                    plan.Execute(bytes, &received_source[0],
                                 &received_parity[0], &restored[0]);
                }, result.selected_parity_null_rejections);
                received_parity[selected] = selected_pointer;
                RequireFilled(restored_storage, 0x3c,
                              result.atomic_rejection_bytes_checked);

                unsigned required_survivor = codec.k;
                for (size_t term_index = plan.offsets.back();
                     term_index > 0; --term_index)
                {
                    const DecodeTerm<Element>& term = plan.terms[term_index - 1];
                    if (!term.parity)
                    {
                        required_survivor = term.index;
                        break;
                    }
                }
                if (required_survivor < codec.k)
                {
                    const void* survivor_pointer =
                        received_source[required_survivor];
                    received_source[required_survivor] = NULL;
                    ExpectRejected([&]() {
                        plan.Execute(bytes, &received_source[0],
                                     &received_parity[0], &restored[0]);
                    }, result.survivor_null_rejections);
                    received_source[required_survivor] = survivor_pointer;
                    RequireFilled(restored_storage, 0x3c,
                                  result.atomic_rejection_bytes_checked);
                }
                if (sizeof(Element) == 2)
                {
                    ExpectRejected([&]() {
                        plan.Execute(3, &received_source[0],
                                     &received_parity[0], &restored[0]);
                    }, result.odd_gf16_rejections);
                    RequireFilled(restored_storage, 0x3c,
                                  result.atomic_rejection_bytes_checked);
                }
            }

            C7TrackedAllocations = 0;
            C7TrackAllocations = true;
            plan.Execute(bytes, &received_source[0], &received_parity[0],
                         &restored[0]);
            C7TrackAllocations = false;
            result.hot_path_allocations += C7TrackedAllocations;
            ++result.decode_executions;
            for (unsigned original = 0; original < k; ++original)
            {
                const bool is_missing = std::binary_search(
                    missing.begin(), missing.end(), original);
                if (restored_storage[original][0] != 0x3c ||
                    restored_storage[original][bytes + 1] != 0x3c)
                    Fail("exact-low decoder changed an output guard");
                if (is_missing)
                {
                    if (memcmp(&restored_storage[original][1],
                               &source[original][1], bytes) != 0)
                        Fail("exact-low decoder restored the wrong original");
                    result.decode_symbol_comparisons += symbols;
                    result.digest = Fnv(result.digest,
                        &restored_storage[original][1], bytes);
                }
                else if (!std::all_of(restored_storage[original].begin() + 1,
                                      restored_storage[original].end() - 1,
                                      [](uint8_t value) { return value == 0x3c; }))
                    Fail("exact-low decoder wrote a surviving original");
            }

            // Re-encoding after recovery must reproduce every parity byte.
            std::vector<const void*> rebuilt_source = source_ptr;
            for (size_t i = 0; i < missing.size(); ++i)
                rebuilt_source[missing[i]] = restored[missing[i]];
            std::vector<std::vector<uint8_t> > rebuilt_parity(
                r, std::vector<uint8_t>(bytes + 2, 0x6d));
            std::vector<void*> rebuilt_ptr(r);
            for (unsigned recovery = 0; recovery < r; ++recovery)
                rebuilt_ptr[recovery] = &rebuilt_parity[recovery][1];
            C7TrackedAllocations = 0;
            C7TrackAllocations = true;
            codec.Encode(bytes, &rebuilt_source[0], &rebuilt_ptr[0]);
            C7TrackAllocations = false;
            result.hot_path_allocations += C7TrackedAllocations;
            ++result.parity_rebuilds;
            for (unsigned recovery = 0; recovery < r; ++recovery)
                if (memcmp(&rebuilt_parity[recovery][1],
                           &parity[recovery][1], bytes) != 0)
                    Fail("re-encoded parity differs after recovery");
        }
    }
}

static Correctness RunCorrectness(const leopard::backend::Ops& backend_ops)
{
    Correctness result = {};
    result.digest = UINT64_C(1469598103934665603);
    const IndependentField<uint8_t> gf8 = IndependentGF8();
    const IndependentField<uint16_t> gf16 = IndependentGF16();
    const unsigned gf8_cases[][2] = {
        { 1, 255 }, { 3, 253 }, { 7, 249 }, { 16, 240 },
        { 64, 192 }, { 127, 129 }, { 192, 64 }, { 248, 8 }, { 255, 1 }
    };
    const unsigned gf16_cases[][2] = {
        { 3, 500 }, { 17, 257 }, { 129, 100 }, { 257, 129 }, { 1000, 17 }
    };
    const size_t gf8_bytes_array[] = { 1, 2, 3, 7, 31, 64, 65, 257 };
    const size_t gf16_bytes_array[] = { 2, 4, 6, 14, 62, 64, 66, 130, 514 };
    const std::vector<size_t> gf8_bytes(
        gf8_bytes_array, gf8_bytes_array + sizeof(gf8_bytes_array) / sizeof(size_t));
    const std::vector<size_t> gf16_bytes(
        gf16_bytes_array,
        gf16_bytes_array + sizeof(gf16_bytes_array) / sizeof(size_t));
    for (unsigned i = 0; i < sizeof(gf8_cases) / sizeof(gf8_cases[0]); ++i)
    {
        ValidateCase<FieldOps<uint8_t> >(
            backend_ops, gf8_cases[i][0], gf8_cases[i][1], i,
            gf8_bytes, gf8, result);
        ++result.gf8_cases;
    }
    for (unsigned i = 0; i < sizeof(gf16_cases) / sizeof(gf16_cases[0]); ++i)
    {
        ValidateCase<FieldOps<uint16_t> >(
            backend_ops, gf16_cases[i][0], gf16_cases[i][1], 100 + i,
            gf16_bytes, gf16, result);
        ++result.gf16_cases;
    }
    ValidateExhaustiveSmallPlans(backend_ops, result);
    ValidateConcurrentBackendContexts(result);
    if (result.hot_path_allocations != 0)
        Fail("exact-low execution allocated in a hot path");
    return result;
}

class AlignedBuffer
{
public:
    AlignedBuffer() : pointer_(NULL), bytes_(0) {}
    ~AlignedBuffer() { free(pointer_); }
    void Reset(size_t bytes)
    {
        free(pointer_);
        pointer_ = NULL;
        bytes_ = 0;
        if (!bytes)
            return;
        if (posix_memalign(&pointer_, leo2_scratch_alignment(), bytes) != 0)
            throw std::bad_alloc();
        bytes_ = bytes;
        memset(pointer_, 0, bytes);
    }
    void* data() { return pointer_; }
    size_t size() const { return bytes_; }
private:
    AlignedBuffer(const AlignedBuffer&);
    AlignedBuffer& operator=(const AlignedBuffer&);
    void* pointer_;
    size_t bytes_;
};

struct Summary
{
    double median;
    double mad;
};

static double Median(std::vector<double> values)
{
    std::sort(values.begin(), values.end());
    const size_t middle = values.size() / 2;
    return values.size() & 1U ? values[middle]
        : (values[middle - 1] + values[middle]) * 0.5;
}

static Summary Summarize(const std::vector<double>& values)
{
    Summary result;
    result.median = Median(values);
    std::vector<double> deviations(values.size());
    for (size_t i = 0; i < values.size(); ++i)
        deviations[i] = std::fabs(values[i] - result.median);
    result.mad = Median(deviations);
    return result;
}

struct BenchmarkSpec
{
    unsigned k;
    unsigned r;
    uint64_t bytes;
    unsigned batch;
    unsigned losses;
    leo2_field exact_field;
    leo2_field padded_field;
};

struct BenchmarkResult
{
    BenchmarkSpec spec;
    Summary exact_setup;
    Summary padded_setup;
    Summary exact_encode;
    Summary padded_encode;
    Summary exact_decode_setup;
    Summary padded_decode_setup;
    Summary exact_decode;
    Summary padded_decode;
    size_t exact_coefficients;
    size_t exact_decode_terms;
    size_t padded_encode_scratch;
    size_t padded_decode_scratch;
    std::vector<double> exact_encode_samples;
    std::vector<double> padded_encode_samples;
    std::vector<double> exact_setup_samples;
    std::vector<double> padded_setup_samples;
    std::vector<double> exact_decode_setup_samples;
    std::vector<double> padded_decode_setup_samples;
    std::vector<double> exact_decode_samples;
    std::vector<double> padded_decode_samples;
};

struct Stripe
{
    std::vector<std::vector<uint8_t> > source;
    std::vector<std::vector<uint8_t> > exact_parity;
    std::vector<std::vector<uint8_t> > padded_parity;
    std::vector<std::vector<uint8_t> > exact_restored;
    std::vector<std::vector<uint8_t> > padded_restored;
    std::vector<const void*> source_ptr;
    std::vector<void*> exact_parity_ptr;
    std::vector<void*> padded_parity_ptr;
    std::vector<const void*> received_source;
    std::vector<const void*> exact_received_parity;
    std::vector<const void*> padded_received_parity;
    std::vector<void*> exact_restored_ptr;
    std::vector<void*> padded_restored_ptr;

    Stripe(unsigned k, unsigned r, size_t bytes, unsigned salt)
        : source(k, std::vector<uint8_t>(bytes))
        , exact_parity(r, std::vector<uint8_t>(bytes))
        , padded_parity(r, std::vector<uint8_t>(bytes))
        , exact_restored(k, std::vector<uint8_t>(bytes))
        , padded_restored(k, std::vector<uint8_t>(bytes))
        , source_ptr(k), exact_parity_ptr(r), padded_parity_ptr(r)
        , received_source(k), exact_received_parity(r)
        , padded_received_parity(r), exact_restored_ptr(k)
        , padded_restored_ptr(k)
    {
        for (unsigned original = 0; original < k; ++original)
        {
            for (size_t byte = 0; byte < bytes; ++byte)
                source[original][byte] = static_cast<uint8_t>(
                    original * 73U + byte * 29U + salt * 11U + 1U);
            source_ptr[original] = &source[original][0];
            received_source[original] = &source[original][0];
            exact_restored_ptr[original] = &exact_restored[original][0];
            padded_restored_ptr[original] = &padded_restored[original][0];
        }
        for (unsigned recovery = 0; recovery < r; ++recovery)
        {
            exact_parity_ptr[recovery] = &exact_parity[recovery][0];
            padded_parity_ptr[recovery] = &padded_parity[recovery][0];
            exact_received_parity[recovery] = &exact_parity[recovery][0];
            padded_received_parity[recovery] = &padded_parity[recovery][0];
        }
    }
};

template <typename Ops>
static BenchmarkResult BenchmarkTyped(
    leo2_context* context, const leopard::backend::Ops& backend_ops,
    const BenchmarkSpec& spec)
{
    typedef ExactCodec<Ops> Codec;
    const unsigned repeats = 7;
    const unsigned warmups = 2;
    BenchmarkResult result = {};
    result.spec = spec;
    result.exact_coefficients = static_cast<size_t>(spec.k) * spec.r;

    for (unsigned sample = 0; sample < repeats; ++sample)
    {
        Clock::time_point begin = Clock::now();
        Codec* exact = new Codec(backend_ops, spec.k, spec.r);
        Clock::time_point end = Clock::now();
        result.exact_setup_samples.push_back(
            std::chrono::duration<double, std::micro>(end - begin).count());
        delete exact;

        leo2_codec* padded = NULL;
        begin = Clock::now();
        Require(leo2_codec_create(context, spec.k, spec.r,
            LEO2_PROFILE_LOW_V1, spec.padded_field, NULL, &padded),
            "C7 padded codec setup");
        end = Clock::now();
        result.padded_setup_samples.push_back(
            std::chrono::duration<double, std::micro>(end - begin).count());
        leo2_codec_destroy(padded);
    }
    result.exact_setup = Summarize(result.exact_setup_samples);
    result.padded_setup = Summarize(result.padded_setup_samples);

    Codec exact(backend_ops, spec.k, spec.r);
    leo2_codec* padded = NULL;
    Require(leo2_codec_create(context, spec.k, spec.r,
        LEO2_PROFILE_LOW_V1, spec.padded_field, NULL, &padded),
        "C7 padded codec create");
    Require(leo2_encode_scratch_size(
        padded, spec.bytes, &result.padded_encode_scratch),
        "C7 padded encode scratch");
    AlignedBuffer encode_scratch;
    encode_scratch.Reset(result.padded_encode_scratch);

    std::vector<Stripe> stripes;
    stripes.reserve(spec.batch);
    for (unsigned stripe = 0; stripe < spec.batch; ++stripe)
        stripes.push_back(Stripe(spec.k, spec.r,
                                 static_cast<size_t>(spec.bytes), stripe));
    for (unsigned stripe = 0; stripe < spec.batch; ++stripe)
    {
        exact.Encode(spec.bytes, &stripes[stripe].source_ptr[0],
                     &stripes[stripe].exact_parity_ptr[0]);
        Require(leo2_encode(padded, spec.bytes,
            &stripes[stripe].source_ptr[0],
            &stripes[stripe].padded_parity_ptr[0], encode_scratch.data(),
            encode_scratch.size()), "C7 padded parity preparation");
    }

    const std::vector<unsigned> missing = MissingSet(spec.k, spec.losses);
    std::vector<uint8_t> parity_present(spec.r, 1);
    ExactDecodePlan<Ops> exact_decode(exact, missing, parity_present);
    result.exact_decode_terms = exact_decode.terms.size();
    std::vector<uint8_t> original_present(spec.k, 1);
    for (size_t i = 0; i < missing.size(); ++i)
        original_present[missing[i]] = 0;
    leo2_decode_plan* padded_decode = NULL;
    Require(leo2_decode_plan_create(padded, &original_present[0],
        &parity_present[0], &padded_decode), "C7 padded decode setup");
    Require(leo2_decode_plan_scratch_size(
        padded_decode, spec.bytes, &result.padded_decode_scratch),
        "C7 padded decode scratch");
    AlignedBuffer decode_scratch;
    decode_scratch.Reset(result.padded_decode_scratch);
    for (unsigned stripe = 0; stripe < spec.batch; ++stripe)
        for (size_t i = 0; i < missing.size(); ++i)
            stripes[stripe].received_source[missing[i]] = NULL;

    for (unsigned sample = 0; sample < repeats; ++sample)
    {
        Clock::time_point begin = Clock::now();
        ExactDecodePlan<Ops>* plan = new ExactDecodePlan<Ops>(
            exact, missing, parity_present);
        Clock::time_point end = Clock::now();
        result.exact_decode_setup_samples.push_back(
            std::chrono::duration<double, std::micro>(end - begin).count());
        delete plan;

        leo2_decode_plan* public_plan = NULL;
        begin = Clock::now();
        Require(leo2_decode_plan_create(padded, &original_present[0],
            &parity_present[0], &public_plan), "C7 padded decode setup sample");
        end = Clock::now();
        result.padded_decode_setup_samples.push_back(
            std::chrono::duration<double, std::micro>(end - begin).count());
        leo2_decode_plan_destroy(public_plan);
    }
    result.exact_decode_setup = Summarize(result.exact_decode_setup_samples);
    result.padded_decode_setup = Summarize(result.padded_decode_setup_samples);

    for (unsigned round = 0; round < warmups + repeats; ++round)
    {
        const bool exact_first = (round & 1U) == 0;
        for (unsigned side = 0; side < 2; ++side)
        {
            const bool run_exact = (side == 0) == exact_first;
            Clock::time_point begin = Clock::now();
            for (unsigned stripe = 0; stripe < spec.batch; ++stripe)
            {
                if (run_exact)
                    exact.Encode(spec.bytes, &stripes[stripe].source_ptr[0],
                                 &stripes[stripe].exact_parity_ptr[0]);
                else
                    Require(leo2_encode(padded, spec.bytes,
                        &stripes[stripe].source_ptr[0],
                        &stripes[stripe].padded_parity_ptr[0],
                        encode_scratch.data(), encode_scratch.size()),
                        "C7 padded encode benchmark");
            }
            Clock::time_point end = Clock::now();
            if (round >= warmups)
            {
                const double micros =
                    std::chrono::duration<double, std::micro>(end - begin).count();
                (run_exact ? result.exact_encode_samples
                           : result.padded_encode_samples).push_back(micros);
            }
        }

        for (unsigned side = 0; side < 2; ++side)
        {
            const bool run_exact = (side == 0) == exact_first;
            Clock::time_point begin = Clock::now();
            for (unsigned stripe = 0; stripe < spec.batch; ++stripe)
            {
                if (run_exact)
                    exact_decode.Execute(spec.bytes,
                        &stripes[stripe].received_source[0],
                        &stripes[stripe].exact_received_parity[0],
                        &stripes[stripe].exact_restored_ptr[0]);
                else
                    Require(leo2_decode_plan_execute(padded_decode, spec.bytes,
                        &stripes[stripe].received_source[0],
                        &stripes[stripe].padded_received_parity[0],
                        &stripes[stripe].padded_restored_ptr[0],
                        decode_scratch.data(), decode_scratch.size()),
                        "C7 padded decode benchmark");
            }
            Clock::time_point end = Clock::now();
            if (round >= warmups)
            {
                const double micros =
                    std::chrono::duration<double, std::micro>(end - begin).count();
                (run_exact ? result.exact_decode_samples
                           : result.padded_decode_samples).push_back(micros);
            }
        }
    }
    result.exact_encode = Summarize(result.exact_encode_samples);
    result.padded_encode = Summarize(result.padded_encode_samples);
    result.exact_decode = Summarize(result.exact_decode_samples);
    result.padded_decode = Summarize(result.padded_decode_samples);
    leo2_decode_plan_destroy(padded_decode);
    leo2_codec_destroy(padded);
    return result;
}

static std::vector<BenchmarkSpec> BenchmarkSpecs(bool smoke)
{
    const BenchmarkSpec all[] = {
        { 3, 253, 64, 8, 3, LEO2_FIELD_GF8, LEO2_FIELD_GF16 },
        { 3, 253, 1024, 8, 3, LEO2_FIELD_GF8, LEO2_FIELD_GF16 },
        { 3, 253, 65536, 1, 3, LEO2_FIELD_GF8, LEO2_FIELD_GF16 },
        { 16, 240, 1024, 8, 8, LEO2_FIELD_GF8, LEO2_FIELD_GF8 },
        { 64, 192, 65536, 1, 8, LEO2_FIELD_GF8, LEO2_FIELD_GF8 },
        { 100, 156, 1024, 8, 8, LEO2_FIELD_GF8, LEO2_FIELD_GF16 },
        { 127, 129, 65536, 1, 8, LEO2_FIELD_GF8, LEO2_FIELD_GF16 },
        { 128, 128, 1024, 8, 8, LEO2_FIELD_GF8, LEO2_FIELD_GF8 },
        { 192, 64, 1024, 8, 8, LEO2_FIELD_GF8, LEO2_FIELD_GF16 },
        { 248, 8, 65536, 1, 8, LEO2_FIELD_GF8, LEO2_FIELD_GF16 },
        { 129, 100, 1024, 8, 8, LEO2_FIELD_GF16, LEO2_FIELD_GF16 },
        { 1000, 200, 1024, 1, 8, LEO2_FIELD_GF16, LEO2_FIELD_GF16 }
    };
    if (smoke)
        return std::vector<BenchmarkSpec>(all, all + 1);
    return std::vector<BenchmarkSpec>(all, all + sizeof(all) / sizeof(all[0]));
}

static void WriteSummary(std::ostream& output, const Summary& summary)
{
    output << "{\"median_us\":" << summary.median
           << ",\"mad_us\":" << summary.mad << '}';
}

static void WriteSamples(std::ostream& output, const std::vector<double>& samples)
{
    output << '[';
    for (size_t i = 0; i < samples.size(); ++i)
    {
        if (i)
            output << ',';
        output << samples[i];
    }
    output << ']';
}

static void WriteJson(
    const std::string& path, const std::string& requested_backend,
    leo2_backend runtime_backend, leo2_backend exact_byte_backend,
    const Correctness& correctness,
    const std::vector<BenchmarkResult>& benchmarks,
    bool production_profile_rejected, const char* timing_scope)
{
    std::ofstream output(path.c_str(), std::ios::binary | std::ios::trunc);
    if (!output)
        Fail("cannot open C7 JSON output");
    const std::vector<unsigned> affinity = ProcessAffinity();
    const char* omp_threads = getenv("OMP_NUM_THREADS");
    const char* omp_dynamic = getenv("OMP_DYNAMIC");
    output << std::setprecision(17);
    output << "{\n"
           << "  \"schema\":\"leopard2-c7-exact-low/v1\",\n"
           << "  \"status\":\"pass\",\n"
           << "  \"profile\":{\"family\":3,\"version\":1,"
              "\"coordinate_map\":1,\"systematic\":\"0..K-1\","
              "\"parity\":\"K..K+R-1\",\"production_enabled\":false},\n"
           << "  \"production_constructor_rejected\":"
           << (production_profile_rejected ? "true" : "false") << ",\n"
           << "  \"timing_scope\":\"" << timing_scope << "\",\n"
           << "  \"requested_backend\":\"" << requested_backend << "\",\n"
           << "  \"runtime_backend\":\"" << BackendName(runtime_backend) << "\",\n"
           << "  \"exact_byte_backend\":\""
           << BackendName(exact_byte_backend) << "\",\n"
           << "  \"affinity\":[";
    for (size_t i = 0; i < affinity.size(); ++i)
    {
        if (i)
            output << ',';
        output << affinity[i];
    }
    output << "],\n"
           << "  \"omp_num_threads\":\""
           << (omp_threads ? omp_threads : "") << "\",\n"
           << "  \"omp_dynamic\":\""
           << (omp_dynamic ? omp_dynamic : "") << "\",\n"
           << "  \"source_sha256\":\"" << LEO2_C7_SOURCE_SHA256 << "\",\n"
           << "  \"core_git_sha\":\"" << LEO2_C7_CORE_GIT_SHA << "\",\n"
           << "  \"library_sha256\":\"" << LEO2_C7_LIBRARY_SHA256 << "\",\n"
           << "  \"sanitizer\":\"" << LEO2_C7_SANITIZER_MODE << "\",\n"
           << "  \"sanitizer_features\":{\"address\":"
           << (LEO2_C7_ASAN_COMPILED ? "true" : "false")
           << ",\"undefined\":"
           << (LEO2_C7_UBSAN_COMPILED ? "true" : "false") << "},\n"
           << "  \"allocation_tracking\":\""
           << LEO2_C7_ALLOCATION_TRACKING_MODE << "\",\n"
           << "  \"correctness\":{\n"
           << "    \"gf8_cases\":" << correctness.gf8_cases << ",\n"
           << "    \"gf16_cases\":" << correctness.gf16_cases << ",\n"
           << "    \"coefficients\":" << correctness.coefficients << ",\n"
           << "    \"gf16_vandermonde_coefficients\":"
           << correctness.gf16_vandermonde_coefficients << ",\n"
           << "    \"encode_executions\":" << correctness.encode_executions << ",\n"
           << "    \"encode_symbol_comparisons\":"
           << correctness.encode_symbol_comparisons << ",\n"
           << "    \"subset_encode_executions\":"
           << correctness.subset_encode_executions << ",\n"
           << "    \"decode_plans\":" << correctness.decode_plans << ",\n"
           << "    \"decode_executions\":" << correctness.decode_executions << ",\n"
           << "    \"decode_symbol_comparisons\":"
           << correctness.decode_symbol_comparisons << ",\n"
           << "    \"maximum_loss_plans\":"
           << correctness.maximum_loss_plans << ",\n"
           << "    \"unavailable_parity_plans\":"
           << correctness.unavailable_parity_plans << ",\n"
           << "    \"no_loss_null_calls\":"
           << correctness.no_loss_null_calls << ",\n"
           << "    \"parity_rebuilds\":" << correctness.parity_rebuilds << ",\n"
           << "    \"odd_gf16_rejections\":"
           << correctness.odd_gf16_rejections << ",\n"
           << "    \"overlap_rejections\":"
           << correctness.overlap_rejections << ",\n"
           << "    \"parity_output_overlap_rejections\":"
           << correctness.parity_output_overlap_rejections << ",\n"
           << "    \"restored_output_overlap_rejections\":"
           << correctness.restored_output_overlap_rejections << ",\n"
           << "    \"restored_input_overlap_rejections\":"
           << correctness.restored_input_overlap_rejections << ",\n"
           << "    \"selected_parity_null_rejections\":"
           << correctness.selected_parity_null_rejections << ",\n"
           << "    \"survivor_null_rejections\":"
           << correctness.survivor_null_rejections << ",\n"
           << "    \"atomic_rejection_bytes_checked\":"
           << correctness.atomic_rejection_bytes_checked << ",\n"
           << "    \"read_only_input_alias_calls\":"
           << correctness.read_only_input_alias_calls << ",\n"
           << "    \"read_only_input_alias_symbol_comparisons\":"
           << correctness.read_only_input_alias_symbol_comparisons << ",\n"
           << "    \"decode_read_only_input_alias_calls\":"
           << correctness.decode_read_only_input_alias_calls << ",\n"
           << "    \"decode_read_only_input_alias_symbol_comparisons\":"
           << correctness.decode_read_only_input_alias_symbol_comparisons
           << ",\n"
           << "    \"detached_plan_executions\":"
           << correctness.detached_plan_executions << ",\n"
           << "    \"detached_plan_symbol_comparisons\":"
           << correctness.detached_plan_symbol_comparisons << ",\n"
           << "    \"concurrent_backend_contexts\":"
           << correctness.concurrent_backend_contexts << ",\n"
           << "    \"concurrent_backend_executions\":"
           << correctness.concurrent_backend_executions << ",\n"
           << "    \"concurrent_backend_trace_calls\":"
           << correctness.concurrent_backend_trace_calls << ",\n"
           << "    \"concurrent_wire_digest_comparisons\":"
           << correctness.concurrent_wire_digest_comparisons << ",\n"
           << "    \"exhaustive_small_plans\":"
           << correctness.exhaustive_small_plans << ",\n"
           << "    \"exhaustive_small_executions\":"
           << correctness.exhaustive_small_executions << ",\n"
           << "    \"exhaustive_small_symbol_comparisons\":"
           << correctness.exhaustive_small_symbol_comparisons << ",\n"
           << "    \"malformed_plan_rejections\":"
           << correctness.malformed_plan_rejections << ",\n"
           << "    \"hot_path_allocations\":"
           << correctness.hot_path_allocations << ",\n"
           << "    \"digest_fnv64\":\"0x" << std::hex
           << correctness.digest << std::dec << "\"\n"
           << "  },\n"
           << "  \"benchmarks\":[\n";
    for (size_t i = 0; i < benchmarks.size(); ++i)
    {
        const BenchmarkResult& cell = benchmarks[i];
        output << "    {\"K\":" << cell.spec.k
               << ",\"R\":" << cell.spec.r
               << ",\"bytes\":" << cell.spec.bytes
               << ",\"batch\":" << cell.spec.batch
               << ",\"losses\":" << cell.spec.losses
               << ",\"exact_field\":" << cell.spec.exact_field
               << ",\"padded_field\":" << cell.spec.padded_field
               << ",\"exact_coefficients\":" << cell.exact_coefficients
               << ",\"exact_decode_terms\":" << cell.exact_decode_terms
               << ",\"padded_encode_scratch\":" << cell.padded_encode_scratch
               << ",\"padded_decode_scratch\":" << cell.padded_decode_scratch
               << ",\"exact_setup\":";
        WriteSummary(output, cell.exact_setup);
        output << ",\"padded_setup\":";
        WriteSummary(output, cell.padded_setup);
        output << ",\"exact_encode\":";
        WriteSummary(output, cell.exact_encode);
        output << ",\"padded_encode\":";
        WriteSummary(output, cell.padded_encode);
        output << ",\"exact_decode_setup\":";
        WriteSummary(output, cell.exact_decode_setup);
        output << ",\"padded_decode_setup\":";
        WriteSummary(output, cell.padded_decode_setup);
        output << ",\"exact_decode\":";
        WriteSummary(output, cell.exact_decode);
        output << ",\"padded_decode\":";
        WriteSummary(output, cell.padded_decode);
        output << ",\"exact_setup_samples_us\":";
        WriteSamples(output, cell.exact_setup_samples);
        output << ",\"padded_setup_samples_us\":";
        WriteSamples(output, cell.padded_setup_samples);
        output << ",\"exact_decode_setup_samples_us\":";
        WriteSamples(output, cell.exact_decode_setup_samples);
        output << ",\"padded_decode_setup_samples_us\":";
        WriteSamples(output, cell.padded_decode_setup_samples);
        output << ",\"exact_encode_samples_us\":";
        WriteSamples(output, cell.exact_encode_samples);
        output << ",\"padded_encode_samples_us\":";
        WriteSamples(output, cell.padded_encode_samples);
        output << ",\"exact_decode_samples_us\":";
        WriteSamples(output, cell.exact_decode_samples);
        output << ",\"padded_decode_samples_us\":";
        WriteSamples(output, cell.padded_decode_samples);
        output << '}' << (i + 1 == benchmarks.size() ? "\n" : ",\n");
    }
    output << "  ]\n}\n";
    if (!output)
        Fail("failed while writing C7 JSON output");
}

} // namespace

int main(int argc, char** argv)
{
    try
    {
        const bool correctness_only = argc == 5 &&
            std::string(argv[4]) == "--correctness-only";
        const bool benchmark_smoke = argc == 5 &&
            std::string(argv[4]) == "--benchmark-smoke";
        if ((argc != 4 && !correctness_only && !benchmark_smoke) ||
            std::string(argv[1]) != "--backend" ||
            (std::string(argv[2]) != "auto" &&
             std::string(argv[2]) != "scalar" &&
             std::string(argv[2]) != "ssse3" &&
             std::string(argv[2]) != "avx2"))
        {
            std::cerr << "usage: " << argv[0]
                      << " --backend NAME OUTPUT.json "
                         "[--correctness-only|--benchmark-smoke]\n";
            return 2;
        }
        leo2_context_options options = {};
        options.struct_size = sizeof(options);
        const std::string requested(argv[2]);
        const leo2_backend requested_backend = ParseBackend(requested);
        options.backend = static_cast<uint32_t>(requested_backend);
        options.thread_count = 1;
        leo2_context* context = NULL;
        Require(leo2_context_create(&options, &context), "C7 context create");
        const leo2_backend runtime = leo2_context_backend(context);
        if (requested != "auto" && requested != BackendName(runtime))
            Fail("runtime backend differs from requested build label");
        leopard::backend::QualificationStatus qualification =
            leopard::backend::QualificationAvailable;
        const leopard::backend::Ops* backend_ops =
            leopard::backend::GetQualifiedOps(
                requested_backend, &qualification);
        if (!backend_ops)
            Fail("exact-low byte backend did not qualify");
        if (runtime != LEO2_BACKEND_NEON && backend_ops->kind != runtime)
            Fail("exact-low byte backend differs from context backend");

        leo2_codec* unsupported = reinterpret_cast<leo2_codec*>(1);
        const leo2_result exact_result = leo2_codec_create(
            context, 3, 253, LEO2_PROFILE_EXACT_EXPERIMENTAL_V1,
            LEO2_FIELD_GF8, NULL, &unsupported);
        const bool production_rejected =
            exact_result == LEO2_UNSUPPORTED && unsupported == NULL;
        if (!production_rejected)
            Fail("production constructor unexpectedly enabled C7 profile");

        const Correctness correctness = RunCorrectness(*backend_ops);
        std::vector<BenchmarkResult> benchmarks;
        if (!correctness_only && std::string(LEO2_C7_SANITIZER_MODE) == "none")
        {
            const std::vector<BenchmarkSpec> specs = BenchmarkSpecs(benchmark_smoke);
            for (size_t i = 0; i < specs.size(); ++i)
            {
                std::cerr << "C7 benchmark " << (i + 1) << '/'
                          << specs.size() << '\n';
                if (specs[i].exact_field == LEO2_FIELD_GF8)
                    benchmarks.push_back(BenchmarkTyped<FieldOps<uint8_t> >(
                        context, *backend_ops, specs[i]));
                else
                    benchmarks.push_back(BenchmarkTyped<FieldOps<uint16_t> >(
                        context, *backend_ops, specs[i]));
            }
        }
        WriteJson(argv[3], requested, runtime, backend_ops->kind,
                  correctness, benchmarks,
                  production_rejected,
                  correctness_only ? "none-correctness-only" :
                  (benchmark_smoke ? "non-authoritative-smoke" :
                   "candidate-authoritative-requires-runner-manifest"));
        leo2_context_destroy(context);
        return 0;
    }
    catch (const std::exception& error)
    {
        std::cerr << "c7_exact_low: " << error.what() << '\n';
        return 1;
    }
}
