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
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <utility>
#include <vector>

#if defined(__linux__)
#include <sched.h>
#endif

#ifndef LEO2_LOCATOR_SOURCE_GIT_SHA
#define LEO2_LOCATOR_SOURCE_GIT_SHA "unknown"
#endif

#ifndef LEO2_LOCATOR_SOURCE_DIRTY
#define LEO2_LOCATOR_SOURCE_DIRTY 1
#endif

namespace {

struct Options
{
    std::string field;
    unsigned n;
    unsigned erasures;
    unsigned calls;
    unsigned iterations;
    unsigned warmup;
    unsigned seed;

    Options()
        : field("gf16")
        , n(1024)
        , erasures(5)
        , calls(16)
        , iterations(101)
        , warmup(9)
        , seed(1)
    {
    }
};

static unsigned ParseUnsigned(const char* text, const char* name)
{
    if (!text[0])
    {
        std::cerr << "invalid " << name << std::endl;
        std::exit(2);
    }
    unsigned value = 0;
    for (const char* next = text; *next; ++next)
    {
        if (*next < '0' || *next > '9')
        {
            std::cerr << "invalid " << name << std::endl;
            std::exit(2);
        }
        const unsigned digit = static_cast<unsigned>(*next - '0');
        if (value > (std::numeric_limits<unsigned>::max() - digit) / 10u)
        {
            std::cerr << "invalid " << name << std::endl;
            std::exit(2);
        }
        value = value * 10u + digit;
    }
    return value;
}

static Options ParseOptions(int argc, char** argv)
{
    Options options;
    for (int i = 1; i < argc; ++i)
    {
        if (!std::strcmp(argv[i], "--help"))
        {
            std::cout
                << "Usage: bench_leopard2_locator [options]\n"
                << "  --field gf8|gf16\n"
                << "  --n POWER_OF_TWO\n"
                << "  --erasures N\n"
                << "  --calls N\n"
                << "  --iterations N\n"
                << "  --warmup N\n"
                << "  --seed N\n";
            std::exit(0);
        }
        if (i + 1 >= argc)
        {
            std::cerr << "missing value for " << argv[i] << std::endl;
            std::exit(2);
        }
        const char* value = argv[++i];
        if (!std::strcmp(argv[i - 1], "--field"))
            options.field = value;
        else if (!std::strcmp(argv[i - 1], "--n"))
            options.n = ParseUnsigned(value, "n");
        else if (!std::strcmp(argv[i - 1], "--erasures"))
            options.erasures = ParseUnsigned(value, "erasures");
        else if (!std::strcmp(argv[i - 1], "--calls"))
            options.calls = ParseUnsigned(value, "calls");
        else if (!std::strcmp(argv[i - 1], "--iterations"))
            options.iterations = ParseUnsigned(value, "iterations");
        else if (!std::strcmp(argv[i - 1], "--warmup"))
            options.warmup = ParseUnsigned(value, "warmup");
        else if (!std::strcmp(argv[i - 1], "--seed"))
            options.seed = ParseUnsigned(value, "seed");
        else
        {
            std::cerr << "unknown option " << argv[i - 1] << std::endl;
            std::exit(2);
        }
    }
    return options;
}

struct Summary
{
    double median_us;
    double mad_us;
    std::vector<double> samples_us;
};

static uint32_t Mix(uint32_t value)
{
    value ^= value >> 16;
    value *= 0x7feb352du;
    value ^= value >> 15;
    value *= 0x846ca68bu;
    return value ^ (value >> 16);
}

template<class Ffe>
static bool Equivalent(const std::vector<Ffe>& a, const std::vector<Ffe>& b)
{
    const Ffe sentinel = std::numeric_limits<Ffe>::max();
    for (size_t i = 0; i < a.size(); ++i)
    {
        const Ffe x = a[i] == sentinel ? 0 : a[i];
        const Ffe y = b[i] == sentinel ? 0 : b[i];
        if (x != y)
            return false;
    }
    return true;
}

template<class Ffe>
static void Warm(
    void (*prepare)(unsigned, const uint8_t*, Ffe*),
    unsigned n,
    const std::vector<uint8_t>& erasures,
    unsigned calls,
    unsigned warmup,
    std::vector<Ffe>& output)
{
    for (unsigned pass = 0; pass < warmup; ++pass)
        for (unsigned call = 0; call < calls; ++call)
            prepare(n, &erasures[0], &output[0]);
}

template<class Ffe>
static double MeasureOne(
    void (*prepare)(unsigned, const uint8_t*, Ffe*),
    unsigned n,
    const std::vector<uint8_t>& erasures,
    unsigned calls,
    std::vector<Ffe>& output)
{
    const std::chrono::steady_clock::time_point start =
        std::chrono::steady_clock::now();
    for (unsigned call = 0; call < calls; ++call)
        prepare(n, &erasures[0], &output[0]);
    const std::chrono::steady_clock::time_point stop =
        std::chrono::steady_clock::now();
    return std::chrono::duration<double, std::micro>(stop - start).count() /
        calls;
}

static Summary Summarize(const std::vector<double>& samples)
{
    std::vector<double> sorted(samples);
    std::sort(sorted.begin(), sorted.end());
    const double median = sorted[sorted.size() / 2];
    std::vector<double> deviations(sorted.size());
    for (size_t i = 0; i < sorted.size(); ++i)
        deviations[i] = sorted[i] > median ?
            sorted[i] - median : median - sorted[i];
    std::sort(deviations.begin(), deviations.end());
    Summary result;
    result.median_us = median;
    result.mad_us = deviations[deviations.size() / 2];
    result.samples_us = samples;
    return result;
}

template<class Ffe>
static std::pair<Summary, Summary> MeasurePair(
    void (*direct)(unsigned, const uint8_t*, Ffe*),
    void (*active)(unsigned, const uint8_t*, Ffe*),
    unsigned n,
    const std::vector<uint8_t>& erasures,
    unsigned calls,
    unsigned iterations,
    unsigned warmup,
    std::vector<Ffe>& direct_output,
    std::vector<Ffe>& active_output)
{
    // Alternate AB then BA.  The resulting global order is ABBA and avoids
    // assigning all thermal/frequency drift to the second implementation.
    for (unsigned pass = 0; pass < warmup; ++pass)
    {
        if ((pass & 1u) == 0)
        {
            Warm(direct, n, erasures, calls, 1, direct_output);
            Warm(active, n, erasures, calls, 1, active_output);
        }
        else
        {
            Warm(active, n, erasures, calls, 1, active_output);
            Warm(direct, n, erasures, calls, 1, direct_output);
        }
    }

    std::vector<double> direct_samples, active_samples;
    direct_samples.reserve(iterations);
    active_samples.reserve(iterations);
    for (unsigned sample = 0; sample < iterations; ++sample)
    {
        if ((sample & 1u) == 0)
        {
            direct_samples.push_back(MeasureOne(
                direct, n, erasures, calls, direct_output));
            active_samples.push_back(MeasureOne(
                active, n, erasures, calls, active_output));
        }
        else
        {
            active_samples.push_back(MeasureOne(
                active, n, erasures, calls, active_output));
            direct_samples.push_back(MeasureOne(
                direct, n, erasures, calls, direct_output));
        }
    }
    return std::make_pair(
        Summarize(direct_samples), Summarize(active_samples));
}

static std::string Environment(const char* name)
{
    const char* value = std::getenv(name);
    return value ? value : "";
}

static const char* CompilerId()
{
#if defined(__clang__)
    return "clang";
#elif defined(__GNUC__)
    return "gcc";
#elif defined(_MSC_VER)
    return "msvc";
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
#elif defined(_MSC_FULL_VER)
#define LEO2_STRINGIFY_INNER(value) #value
#define LEO2_STRINGIFY(value) LEO2_STRINGIFY_INNER(value)
    return LEO2_STRINGIFY(_MSC_FULL_VER);
#else
    return "unknown";
#endif
}

static std::vector<unsigned> AllowedCpus()
{
    std::vector<unsigned> result;
#if defined(__linux__)
    cpu_set_t set;
    CPU_ZERO(&set);
    if (sched_getaffinity(0, sizeof(set), &set) == 0)
    {
        for (unsigned cpu = 0; cpu < CPU_SETSIZE; ++cpu)
            if (CPU_ISSET(cpu, &set))
                result.push_back(cpu);
    }
#endif
    return result;
}

static void WriteJsonString(std::ostream& output, const std::string& value)
{
    output << '"';
    for (size_t i = 0; i < value.size(); ++i)
    {
        const unsigned char c = static_cast<unsigned char>(value[i]);
        if (c == '"' || c == '\\')
            output << '\\' << c;
        else if (c == '\n')
            output << "\\n";
        else if (c == '\r')
            output << "\\r";
        else if (c == '\t')
            output << "\\t";
        else if (c >= 0x20)
            output << c;
        else
            output << "\\u00" << std::hex << std::setw(2)
                   << std::setfill('0') << static_cast<unsigned>(c)
                   << std::dec << std::setfill(' ');
    }
    output << '"';
}

static void WriteSamples(const std::vector<double>& samples)
{
    std::cout << '[';
    for (size_t i = 0; i < samples.size(); ++i)
    {
        if (i)
            std::cout << ',';
        std::cout << samples[i];
    }
    std::cout << ']';
}

template<class Ffe>
static int Run(
    const Options& options,
    unsigned order,
    bool (*initialize)(),
    void (*direct)(unsigned, const uint8_t*, Ffe*),
    void (*active)(unsigned, const uint8_t*, Ffe*),
    void (*reference)(unsigned, const uint8_t*, Ffe*),
    bool (*direct_preferred)(unsigned, unsigned))
{
    if (options.n < 2 || options.n > order ||
        (options.n & (options.n - 1)) != 0 ||
        options.erasures > options.n || options.calls == 0 ||
        options.iterations == 0)
    {
        std::cerr << "invalid benchmark dimensions" << std::endl;
        return 2;
    }
    if (!initialize())
    {
        std::cerr << "field initialization failed" << std::endl;
        return 1;
    }

    std::vector<uint8_t> erasures(options.n, 0);
    const uint32_t multiplier = Mix(options.seed) | 1u;
    const uint32_t offset = Mix(options.seed ^ 0x6d2b79f5u);
    for (unsigned i = 0; i < options.erasures; ++i)
        erasures[(i * multiplier + offset) & (options.n - 1)] = 1;

    std::vector<Ffe> direct_output(options.n), active_output(options.n),
        expected(options.n);
    direct(options.n, &erasures[0], &direct_output[0]);
    active(options.n, &erasures[0], &active_output[0]);
    reference(options.n, &erasures[0], &expected[0]);
    if (!Equivalent(direct_output, expected) ||
        !Equivalent(active_output, expected))
    {
        std::cerr << "locator oracle mismatch" << std::endl;
        return 1;
    }

    const std::pair<Summary, Summary> measured = MeasurePair(
        direct, active, options.n, erasures, options.calls,
        options.iterations, options.warmup, direct_output, active_output);
    const Summary& direct_time = measured.first;
    const Summary& active_time = measured.second;
    const std::vector<unsigned> allowed_cpus = AllowedCpus();

    std::cout << std::setprecision(17)
              << "{\n"
              << "  \"schema\": \"leopard2-locator-benchmark-v2\",\n"
              << "  \"build\": { \"source_git_sha\": ";
    WriteJsonString(std::cout, LEO2_LOCATOR_SOURCE_GIT_SHA);
    std::cout << ", \"source_dirty_at_configure\": "
              << (LEO2_LOCATOR_SOURCE_DIRTY ? "true" : "false")
              << ", \"compiler\": ";
    WriteJsonString(std::cout, CompilerId());
    std::cout << ", \"compiler_version\": ";
    WriteJsonString(std::cout, CompilerVersion());
    std::cout << " },\n"
              << "  \"runtime\": { \"allowed_cpus\": [";
    for (size_t i = 0; i < allowed_cpus.size(); ++i)
    {
        if (i)
            std::cout << ',';
        std::cout << allowed_cpus[i];
    }
    std::cout << "], \"openmp_macro\": ";
#if defined(_OPENMP)
    std::cout << _OPENMP;
#else
    std::cout << 0;
#endif
    std::cout << ", \"omp_num_threads\": ";
    WriteJsonString(std::cout, Environment("OMP_NUM_THREADS"));
    std::cout << ", \"omp_dynamic\": ";
    WriteJsonString(std::cout, Environment("OMP_DYNAMIC"));
    std::cout << " },\n"
              << "  \"field\": \"" << options.field << "\",\n"
              << "  \"parent_n\": " << options.n << ",\n"
              << "  \"erasure_count\": " << options.erasures << ",\n"
              << "  \"calls_per_sample\": " << options.calls << ",\n"
              << "  \"iterations\": " << options.iterations << ",\n"
              << "  \"warmup_calls\": "
              << static_cast<uint64_t>(options.warmup) * options.calls
              << ",\n"
              << "  \"seed\": " << options.seed << ",\n"
              << "  \"measurement_order\": \"ABBA\",\n"
              << "  \"dispatcher\": \""
              << (direct_preferred(options.n, options.erasures) ?
                    "direct" : "active_walsh") << "\",\n"
              << "  \"direct\": { \"median_us\": "
              << direct_time.median_us << ", \"mad_us\": "
              << direct_time.mad_us << ", \"samples_us\": ";
    WriteSamples(direct_time.samples_us);
    std::cout << " },\n"
              << "  \"active_walsh\": { \"median_us\": "
              << active_time.median_us << ", \"mad_us\": "
              << active_time.mad_us << ", \"samples_us\": ";
    WriteSamples(active_time.samples_us);
    std::cout << " }\n"
              << "}\n";
    return 0;
}

} // namespace

int main(int argc, char** argv)
{
    const Options options = ParseOptions(argc, argv);
    if (options.field == "gf8")
    {
#ifdef LEO_HAS_FF8
        return Run<leopard::ff8::ffe_t>(options, leopard::ff8::kOrder,
            leopard::ff8::Initialize, leopard::ff8::PrepareDecodeDirect,
            leopard::ff8::PrepareDecodeWalshActive,
            leopard::ff8::PrepareDecodeWalshReference,
            leopard::ff8::IsDirectLocatorPreferred);
#else
        std::cerr << "GF8 is unavailable in this build" << std::endl;
        return 2;
#endif
    }
    if (options.field == "gf16")
    {
#ifdef LEO_HAS_FF16
        return Run<leopard::ff16::ffe_t>(options, leopard::ff16::kOrder,
            leopard::ff16::Initialize, leopard::ff16::PrepareDecodeDirect,
            leopard::ff16::PrepareDecodeWalshActive,
            leopard::ff16::PrepareDecodeWalshReference,
            leopard::ff16::IsDirectLocatorPreferred);
#else
        std::cerr << "GF16 is unavailable in this build" << std::endl;
        return 2;
#endif
    }
    std::cerr << "field must be gf8 or gf16" << std::endl;
    return 2;
}
