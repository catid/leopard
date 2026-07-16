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

// Default-off C6 experiment.  This file is not part of the root CMake graph.

#include "LeopardFF8.h"
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

#ifndef LEO2_C6_SOURCE_SHA256
#define LEO2_C6_SOURCE_SHA256 "unbound"
#endif
#ifndef LEO2_C6_CORE_GIT_SHA
#define LEO2_C6_CORE_GIT_SHA "unbound"
#endif
#ifndef LEO2_C6_LIBRARY_SHA256
#define LEO2_C6_LIBRARY_SHA256 "unbound"
#endif
#ifndef LEO2_C6_SANITIZER_MODE
#define LEO2_C6_SANITIZER_MODE "none"
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

static unsigned Parent(leo2_profile profile, unsigned k, unsigned r)
{
    return profile == LEO2_PROFILE_LEGACY_HIGH_V1
        ? NextPow2(k + NextPow2(r))
        : NextPow2(NextPow2(k) + r);
}

static const char* ProfileName(leo2_profile profile)
{
    return profile == LEO2_PROFILE_LEGACY_HIGH_V1
        ? "legacy_high_v1" : "low_v1";
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
        if (posix_memalign(&pointer_, 64, bytes) != 0)
            throw std::bad_alloc();
        bytes_ = bytes;
        memset(pointer_, 0, bytes);
    }
    uint8_t* data() { return static_cast<uint8_t*>(pointer_); }
    const uint8_t* data() const { return static_cast<const uint8_t*>(pointer_); }
    size_t size() const { return bytes_; }
private:
    AlignedBuffer(const AlignedBuffer&);
    AlignedBuffer& operator=(const AlignedBuffer&);
    void* pointer_;
    size_t bytes_;
};

struct IndependentGF8
{
    uint8_t to_polynomial[256];
    uint8_t to_coordinate[256];

    IndependentGF8()
    {
        static const unsigned basis[8] = { 1, 214, 152, 146, 86, 200, 88, 230 };
        memset(to_coordinate, 0xff, sizeof(to_coordinate));
        for (unsigned coordinate = 0; coordinate < 256; ++coordinate)
        {
            unsigned polynomial = 0;
            for (unsigned bit = 0; bit < 8; ++bit)
                if (coordinate & (1U << bit))
                    polynomial ^= basis[bit];
            if (polynomial >= 256 || to_coordinate[polynomial] != 0xff)
                Fail("legacy GF8 coordinate basis is not independent");
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
            if (b & 1U)
                product ^= a;
            b >>= 1;
            a <<= 1;
        }
        for (int bit = 14; bit >= 8; --bit)
            if (product & (1U << bit))
                product ^= 0x11dU << (bit - 8);
        return to_coordinate[product];
    }

    uint8_t Power(uint8_t value, unsigned exponent) const
    {
        uint8_t result = 1;
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

    uint8_t Inverse(uint8_t value) const
    {
        if (!value)
            Fail("independent inverse of zero");
        return Power(value, 254);
    }
};

struct ExactPlan
{
    unsigned k;
    unsigned r;
    std::vector<uint8_t> logs;

    ExactPlan(unsigned original_count, unsigned recovery_count)
        : k(original_count), r(recovery_count), logs(k * r)
    {
        std::vector<uint8_t> weights(k);
        for (unsigned i = 0; i < k; ++i)
        {
            uint8_t denominator = 1;
            for (unsigned other = 0; other < k; ++other)
                if (other != i)
                    denominator = leopard::ff8::MultiplyElements(
                        denominator, static_cast<uint8_t>(i ^ other));
            weights[i] = leopard::ff8::InverseElement(denominator);
        }
        for (unsigned parity = 0; parity < r; ++parity)
        {
            const uint8_t point = static_cast<uint8_t>(k + parity);
            uint8_t vanishing = 1;
            for (unsigned i = 0; i < k; ++i)
                vanishing = leopard::ff8::MultiplyElements(
                    vanishing, static_cast<uint8_t>(point ^ i));
            for (unsigned i = 0; i < k; ++i)
            {
                uint8_t coefficient = leopard::ff8::MultiplyElements(
                    vanishing,
                    leopard::ff8::InverseElement(static_cast<uint8_t>(point ^ i)));
                coefficient = leopard::ff8::MultiplyElements(coefficient, weights[i]);
                if (!coefficient)
                    Fail("exact GF8 generator unexpectedly contains zero");
                logs[parity * k + i] = leopard::ff8::ElementLog(coefficient);
            }
        }
    }

    void Encode(uint64_t bytes, const void* const* original, void* const* recovery) const
    {
        for (unsigned parity = 0; parity < r; ++parity)
        {
            uint8_t* output = static_cast<uint8_t*>(recovery[parity]);
            const uint8_t first = logs[parity * k];
            if (first == 0)
                memcpy(output, original[0], static_cast<size_t>(bytes));
            else
                leopard::ff8::MultiplyBytes(output, original[0], first, bytes);
            for (unsigned i = 1; i < k; ++i)
            {
                leopard::ff8::MultiplyAddBytes(
                    output, original[i], logs[parity * k + i], bytes);
            }
        }
    }
};

static std::vector<uint8_t> IndependentRows(
    const IndependentGF8& field, unsigned k, unsigned r)
{
    std::vector<uint8_t> rows(k * r);
    std::vector<uint8_t> weights(k);
    for (unsigned i = 0; i < k; ++i)
    {
        uint8_t denominator = 1;
        for (unsigned other = 0; other < k; ++other)
            if (other != i)
                denominator = field.Multiply(
                    denominator, static_cast<uint8_t>(i ^ other));
        weights[i] = field.Inverse(denominator);
    }
    for (unsigned parity = 0; parity < r; ++parity)
    {
        const uint8_t point = static_cast<uint8_t>(k + parity);
        uint8_t vanishing = 1;
        for (unsigned i = 0; i < k; ++i)
            vanishing = field.Multiply(vanishing, static_cast<uint8_t>(point ^ i));
        for (unsigned i = 0; i < k; ++i)
        {
            rows[parity * k + i] = field.Multiply(
                field.Multiply(vanishing,
                    field.Inverse(static_cast<uint8_t>(point ^ i))),
                weights[i]);
        }
    }
    return rows;
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

struct Correctness
{
    unsigned cases;
    uint64_t coefficients;
    uint64_t byte_comparisons;
    uint64_t digest;
};

static Correctness RunCorrectness()
{
    static const unsigned cases[][2] = {
        { 1, 129 }, { 3, 129 }, { 17, 129 }, { 65, 129 },
        { 127, 129 }, { 129, 1 }, { 129, 3 }, { 129, 17 },
        { 129, 65 }, { 129, 127 }, { 31, 193 }, { 193, 31 }
    };
    IndependentGF8 field;
    Correctness result = { 0, 0, 0, UINT64_C(1469598103934665603) };
    const size_t bytes = 257;
    for (unsigned case_i = 0; case_i < sizeof(cases) / sizeof(cases[0]); ++case_i)
    {
        const unsigned k = cases[case_i][0];
        const unsigned r = cases[case_i][1];
        ExactPlan plan(k, r);
        const std::vector<uint8_t> rows = IndependentRows(field, k, r);
        for (unsigned parity = 0; parity < r; ++parity)
            for (unsigned i = 0; i < k; ++i)
            {
                const uint8_t coefficient = leopard::ff8::MultiplyLogElement(
                    1, plan.logs[parity * k + i]);
                if (coefficient != rows[parity * k + i])
                    Fail("production and independent exact rows disagree");
            }

        std::vector<std::vector<uint8_t> > source(k, std::vector<uint8_t>(bytes + 2));
        std::vector<std::vector<uint8_t> > output(r, std::vector<uint8_t>(bytes + 2, 0xa5));
        std::vector<const void*> source_ptr(k);
        std::vector<void*> output_ptr(r);
        for (unsigned i = 0; i < k; ++i)
        {
            for (size_t byte = 0; byte < bytes; ++byte)
                source[i][byte + 1] = static_cast<uint8_t>(
                    i * 73U + byte * 29U + case_i * 11U + 1U);
            source_ptr[i] = &source[i][1];
        }
        for (unsigned parity = 0; parity < r; ++parity)
            output_ptr[parity] = &output[parity][1];
        plan.Encode(bytes, &source_ptr[0], &output_ptr[0]);
        for (unsigned parity = 0; parity < r; ++parity)
        {
            if (output[parity][0] != 0xa5 || output[parity][bytes + 1] != 0xa5)
                Fail("exact encoder changed guard byte");
            for (size_t byte = 0; byte < bytes; ++byte)
            {
                uint8_t expected = 0;
                for (unsigned i = 0; i < k; ++i)
                    expected ^= field.Multiply(
                        rows[parity * k + i], source[i][byte + 1]);
                if (output[parity][byte + 1] != expected)
                    Fail("exact encoder differs from independent byte algebra");
            }
            result.digest = Fnv(result.digest, &output[parity][1], bytes);
        }
        ++result.cases;
        result.coefficients += static_cast<uint64_t>(k) * r;
        result.byte_comparisons += static_cast<uint64_t>(r) * bytes;
    }
    return result;
}

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

struct Geometry
{
    leo2_profile profile;
    unsigned k;
    unsigned r;
};

struct CellSpec
{
    Geometry geometry;
    uint64_t bytes;
    unsigned batch;
    unsigned reuse;
};

struct Stripe
{
    AlignedBuffer source_storage;
    AlignedBuffer exact_storage;
    AlignedBuffer padded_storage;
    AlignedBuffer scratch;
    std::vector<const void*> source;
    std::vector<void*> exact;
    std::vector<void*> padded;
    leo2_encode_batch_item item;
};

struct CellResult
{
    CellSpec spec;
    Summary exact_setup;
    Summary padded_setup;
    Summary exact_execution;
    Summary padded_execution;
    size_t padded_scratch;
    uint64_t exact_output_digest;
    uint64_t padded_output_digest;
};

static leo2_codec* CreateCodec(leo2_context* context, const Geometry& geometry)
{
    leo2_codec_options options = {};
    options.struct_size = sizeof(options);
    leo2_codec* codec = NULL;
    Require(leo2_codec_create(context, geometry.k, geometry.r,
        geometry.profile, LEO2_FIELD_GF16, &options, &codec), "codec create");
    if (leo2_codec_parent_count(codec) != Parent(
            geometry.profile, geometry.k, geometry.r) ||
        leo2_codec_field(codec) != LEO2_FIELD_GF16)
    {
        leo2_codec_destroy(codec);
        Fail("padded baseline geometry differs from request");
    }
    return codec;
}

static std::vector<CellSpec> BenchmarkSpecs()
{
    static const Geometry geometries[] = {
        { LEO2_PROFILE_LEGACY_HIGH_V1, 2, 129 },
        { LEO2_PROFILE_LEGACY_HIGH_V1, 3, 129 },
        { LEO2_PROFILE_LEGACY_HIGH_V1, 3, 130 },
        { LEO2_PROFILE_LEGACY_HIGH_V1, 3, 131 },
        { LEO2_PROFILE_LEGACY_HIGH_V1, 4, 129 },
        { LEO2_PROFILE_LEGACY_HIGH_V1, 17, 129 },
        { LEO2_PROFILE_LEGACY_HIGH_V1, 17, 130 },
        { LEO2_PROFILE_LOW_V1, 129, 2 },
        { LEO2_PROFILE_LOW_V1, 129, 3 },
        { LEO2_PROFILE_LOW_V1, 130, 3 },
        { LEO2_PROFILE_LOW_V1, 131, 3 },
        { LEO2_PROFILE_LOW_V1, 129, 4 },
        { LEO2_PROFILE_LOW_V1, 129, 17 },
        { LEO2_PROFILE_LOW_V1, 130, 17 }
    };
    std::vector<CellSpec> result;
    for (unsigned i = 0; i < sizeof(geometries) / sizeof(geometries[0]); ++i)
    {
        result.push_back(CellSpec{ geometries[i], 64, 8, 8 });
        result.push_back(CellSpec{ geometries[i], 1024, 8, 4 });
        result.push_back(CellSpec{ geometries[i], 65536, 1, 1 });
        if (geometries[i].k <= 3 || geometries[i].r <= 3)
            result.push_back(CellSpec{ geometries[i], 1048576, 1, 1 });
    }
    return result;
}

static void FillStripe(Stripe& stripe, const CellSpec& spec,
                       size_t scratch_bytes, unsigned stripe_index)
{
    const size_t source_bytes = static_cast<size_t>(spec.geometry.k) * spec.bytes;
    const size_t output_bytes = static_cast<size_t>(spec.geometry.r) * spec.bytes;
    stripe.source_storage.Reset(source_bytes);
    stripe.exact_storage.Reset(output_bytes);
    stripe.padded_storage.Reset(output_bytes);
    stripe.scratch.Reset(scratch_bytes);
    stripe.source.resize(spec.geometry.k);
    stripe.exact.resize(spec.geometry.r);
    stripe.padded.resize(spec.geometry.r);
    for (size_t i = 0; i < source_bytes; ++i)
        stripe.source_storage.data()[i] = static_cast<uint8_t>(
            i * 73U + stripe_index * 41U + spec.geometry.k * 11U + spec.geometry.r);
    for (unsigned i = 0; i < spec.geometry.k; ++i)
        stripe.source[i] = stripe.source_storage.data() + static_cast<size_t>(i) * spec.bytes;
    for (unsigned i = 0; i < spec.geometry.r; ++i)
    {
        stripe.exact[i] = stripe.exact_storage.data() + static_cast<size_t>(i) * spec.bytes;
        stripe.padded[i] = stripe.padded_storage.data() + static_cast<size_t>(i) * spec.bytes;
    }
    stripe.item.shard_bytes = spec.bytes;
    stripe.item.original = &stripe.source[0];
    stripe.item.recovery = &stripe.padded[0];
    stripe.item.scratch = stripe.scratch.data();
    stripe.item.scratch_bytes = scratch_bytes;
}

static CellResult BenchmarkCell(leo2_context* context, const CellSpec& spec)
{
    CellResult result;
    result.spec = spec;
    result.exact_output_digest = UINT64_C(1469598103934665603);
    result.padded_output_digest = UINT64_C(1469598103934665603);

    std::vector<double> exact_setup, padded_setup;
    for (unsigned sample = 0; sample < 9; ++sample)
    {
        Clock::time_point begin = Clock::now();
        ExactPlan* exact = new ExactPlan(spec.geometry.k, spec.geometry.r);
        Clock::time_point end = Clock::now();
        exact_setup.push_back(std::chrono::duration<double, std::micro>(end - begin).count());
        delete exact;

        begin = Clock::now();
        leo2_codec* codec = CreateCodec(context, spec.geometry);
        end = Clock::now();
        padded_setup.push_back(std::chrono::duration<double, std::micro>(end - begin).count());
        leo2_codec_destroy(codec);
    }
    result.exact_setup = Summarize(exact_setup);
    result.padded_setup = Summarize(padded_setup);

    ExactPlan exact(spec.geometry.k, spec.geometry.r);
    leo2_codec* codec = CreateCodec(context, spec.geometry);
    size_t scratch_bytes = 0;
    Require(leo2_encode_scratch_size(codec, spec.bytes, &scratch_bytes),
            "scratch query");
    result.padded_scratch = scratch_bytes;

    std::vector<std::unique_ptr<Stripe> > stripes;
    std::vector<leo2_encode_batch_item> items(spec.batch);
    for (unsigned stripe_i = 0; stripe_i < spec.batch; ++stripe_i)
    {
        stripes.push_back(std::unique_ptr<Stripe>(new Stripe));
        FillStripe(*stripes.back(), spec, scratch_bytes, stripe_i);
        items[stripe_i] = stripes.back()->item;
    }

    for (unsigned warmup = 0; warmup < 2; ++warmup)
    {
        for (unsigned stripe_i = 0; stripe_i < spec.batch; ++stripe_i)
            exact.Encode(spec.bytes, &stripes[stripe_i]->source[0],
                         &stripes[stripe_i]->exact[0]);
        Require(leo2_encode_batch(codec, &items[0], items.size()), "padded warmup");
    }

    std::vector<double> exact_samples, padded_samples;
    for (unsigned sample = 0; sample < 9; ++sample)
    {
        double measured[2] = { 0, 0 };
        for (unsigned order = 0; order < 2; ++order)
        {
            const bool run_exact = ((sample + order) & 1U) == 0;
            const Clock::time_point begin = Clock::now();
            for (unsigned call = 0; call < spec.reuse; ++call)
            {
                if (run_exact)
                {
                    for (unsigned stripe_i = 0; stripe_i < spec.batch; ++stripe_i)
                        exact.Encode(spec.bytes, &stripes[stripe_i]->source[0],
                                     &stripes[stripe_i]->exact[0]);
                }
                else
                {
                    Require(leo2_encode_batch(codec, &items[0], items.size()),
                            "padded benchmark");
                }
            }
            const Clock::time_point end = Clock::now();
            measured[run_exact ? 0 : 1] =
                std::chrono::duration<double, std::micro>(end - begin).count() /
                spec.reuse;
        }
        exact_samples.push_back(measured[0]);
        padded_samples.push_back(measured[1]);
    }
    result.exact_execution = Summarize(exact_samples);
    result.padded_execution = Summarize(padded_samples);

    for (unsigned stripe_i = 0; stripe_i < spec.batch; ++stripe_i)
    {
        result.exact_output_digest = Fnv(result.exact_output_digest,
            stripes[stripe_i]->exact_storage.data(), stripes[stripe_i]->exact_storage.size());
        result.padded_output_digest = Fnv(result.padded_output_digest,
            stripes[stripe_i]->padded_storage.data(), stripes[stripe_i]->padded_storage.size());
    }
    leo2_codec_destroy(codec);
    return result;
}

static std::vector<int> Affinity()
{
    std::vector<int> output;
#if defined(__linux__)
    cpu_set_t set;
    CPU_ZERO(&set);
    if (sched_getaffinity(0, sizeof(set), &set) == 0)
        for (int cpu = 0; cpu < CPU_SETSIZE; ++cpu)
            if (CPU_ISSET(cpu, &set)) output.push_back(cpu);
#endif
    return output;
}

static std::string Hex64(uint64_t value)
{
    std::ostringstream stream;
    stream << "0x" << std::hex << std::setw(16) << std::setfill('0') << value;
    return stream.str();
}

static std::string CpuModel()
{
    std::ifstream input("/proc/cpuinfo");
    std::string line;
    while (std::getline(input, line))
    {
        const std::string key = "model name";
        if (line.compare(0, key.size(), key) != 0)
            continue;
        const size_t colon = line.find(':');
        if (colon != std::string::npos)
        {
            size_t begin = colon + 1;
            while (begin < line.size() && line[begin] == ' ') ++begin;
            return line.substr(begin);
        }
    }
    return "unavailable";
}

static const char* CompilerName()
{
#if defined(__clang__)
    return "clang";
#elif defined(__GNUC__)
    return "gcc";
#else
    return "unknown";
#endif
}

static const char* CompilerVersion()
{
#if defined(__clang__)
    return __clang_version__;
#elif defined(__GNUC__)
    return __VERSION__;
#else
    return "unknown";
#endif
}

static void WriteSummary(std::ostream& out, const Summary& value)
{
    out << "{\"median_us\":" << std::fixed << std::setprecision(6) << value.median
        << ",\"mad_us\":" << value.mad << "}";
}

static void WriteJson(const std::string& path, const char* requested_backend,
                      leo2_backend runtime_backend, const Correctness& correctness,
                      const std::vector<CellResult>& cells)
{
    std::ofstream file(path.c_str(), std::ios::binary | std::ios::trunc);
    if (!file)
        Fail("cannot open JSON output");
    const std::vector<int> affinity = Affinity();
    const char* omp_env = getenv("OMP_NUM_THREADS");
    file << "{\n  \"schema\":\"leopard2-c6-cpp/v1\",\n"
         << "  \"source_sha256\":\"" << LEO2_C6_SOURCE_SHA256 << "\",\n"
         << "  \"core_git_sha\":\"" << LEO2_C6_CORE_GIT_SHA << "\",\n"
         << "  \"library_sha256\":\"" << LEO2_C6_LIBRARY_SHA256 << "\",\n"
         << "  \"sanitizer_mode\":\"" << LEO2_C6_SANITIZER_MODE << "\",\n"
         << "  \"compiler\":\"" << CompilerName() << "\",\n"
         << "  \"compiler_version\":\"" << CompilerVersion() << "\",\n"
         << "  \"cpu_model\":\"" << CpuModel() << "\",\n"
         << "  \"requested_backend\":\"" << requested_backend << "\",\n"
         << "  \"runtime_backend\":\"" << BackendName(runtime_backend) << "\",\n"
         << "  \"affinity\":[";
    for (size_t i = 0; i < affinity.size(); ++i)
    {
        if (i) file << ',';
        file << affinity[i];
    }
    file << "],\n  \"omp_num_threads\":\"" << (omp_env ? omp_env : "") << "\",\n"
         << "  \"openmp_max_threads\":";
#if defined(_OPENMP)
    file << omp_get_max_threads();
#else
    file << 1;
#endif
    file << ",\n  \"correctness\":{\"cases\":" << correctness.cases
         << ",\"coefficients\":" << correctness.coefficients
         << ",\"byte_comparisons\":" << correctness.byte_comparisons
         << ",\"digest\":\"" << Hex64(correctness.digest) << "\"},\n"
         << "  \"cells\":[\n";
    for (size_t index = 0; index < cells.size(); ++index)
    {
        const CellResult& cell = cells[index];
        const uint64_t terms = static_cast<uint64_t>(cell.spec.geometry.k) *
            cell.spec.geometry.r * cell.spec.batch;
        const uint64_t payload = terms * cell.spec.bytes;
        const double ratio = cell.padded_execution.median / cell.exact_execution.median;
        const double credible = 100.0 * ((cell.padded_execution.median -
            cell.padded_execution.mad) / (cell.exact_execution.median +
            cell.exact_execution.mad) - 1.0);
        file << "    {\"profile\":\"" << ProfileName(cell.spec.geometry.profile)
             << "\",\"K\":" << cell.spec.geometry.k
             << ",\"R\":" << cell.spec.geometry.r
             << ",\"parent\":" << Parent(cell.spec.geometry.profile,
                                              cell.spec.geometry.k,
                                              cell.spec.geometry.r)
             << ",\"bytes\":" << cell.spec.bytes
             << ",\"batch\":" << cell.spec.batch
             << ",\"reuse\":" << cell.spec.reuse
             << ",\"exact_setup\":";
        WriteSummary(file, cell.exact_setup);
        file << ",\"padded_setup\":";
        WriteSummary(file, cell.padded_setup);
        file << ",\"exact_execution\":";
        WriteSummary(file, cell.exact_execution);
        file << ",\"padded_execution\":";
        WriteSummary(file, cell.padded_execution);
        file << ",\"padded_over_exact_ratio\":" << std::fixed << std::setprecision(6)
             << ratio << ",\"credible_gain_percent\":" << credible
             << ",\"exact_table_bytes\":"
             << static_cast<uint64_t>(cell.spec.geometry.k) * cell.spec.geometry.r
             << ",\"exact_execution_scratch_bytes\":0"
             << ",\"padded_execution_scratch_bytes\":" << cell.padded_scratch
             << ",\"exact_fixed_multiply_terms\":" << terms
             << ",\"exact_coefficient_payload_bytes\":" << payload
             << ",\"input_bytes\":"
             << static_cast<uint64_t>(cell.spec.geometry.k) * cell.spec.bytes * cell.spec.batch
             << ",\"output_bytes\":"
             << static_cast<uint64_t>(cell.spec.geometry.r) * cell.spec.bytes * cell.spec.batch
             << ",\"exact_output_digest\":\"" << Hex64(cell.exact_output_digest)
             << "\",\"padded_output_digest\":\"" << Hex64(cell.padded_output_digest)
             << "\"}" << (index + 1 == cells.size() ? "\n" : ",\n");
    }
    file << "  ]\n}\n";
    if (!file)
        Fail("failed while writing JSON output");
}

} // namespace

int main(int argc, char** argv)
{
    try
    {
        const bool correctness_only = argc == 5 &&
            std::string(argv[4]) == "--correctness-only";
        if ((argc != 4 && !correctness_only) ||
            std::string(argv[1]) != "--backend" ||
            (std::string(argv[2]) != "auto" && std::string(argv[2]) != "scalar" &&
             std::string(argv[2]) != "ssse3" && std::string(argv[2]) != "avx2"))
        {
            std::cerr << "usage: " << argv[0]
                      << " --backend NAME OUTPUT.json [--correctness-only]\n";
            return 2;
        }
        leo2_context_options context_options = {};
        context_options.struct_size = sizeof(context_options);
        context_options.thread_count = 1;
        leo2_context* context = NULL;
        Require(leo2_context_create(&context_options, &context), "context create");
        const leo2_backend runtime = leo2_context_backend(context);
        const std::string requested(argv[2]);
        if (requested != "auto" && requested != BackendName(runtime))
            Fail("runtime backend differs from requested build label");

        const Correctness correctness = RunCorrectness();
        std::vector<CellResult> cells;
        if (!correctness_only && std::string(LEO2_C6_SANITIZER_MODE) == "none")
        {
            const std::vector<CellSpec> specs = BenchmarkSpecs();
            for (size_t i = 0; i < specs.size(); ++i)
            {
                std::cerr << "C6 benchmark " << (i + 1) << '/' << specs.size() << '\n';
                cells.push_back(BenchmarkCell(context, specs[i]));
            }
        }
        WriteJson(argv[3], argv[2], runtime, correctness, cells);
        leo2_context_destroy(context);
        return 0;
    }
    catch (const std::exception& error)
    {
        std::cerr << "c6_gf256: " << error.what() << '\n';
        return 1;
    }
}
