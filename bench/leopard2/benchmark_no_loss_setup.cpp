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

#include "Leopard2Direct.h"
#include "leopard2.h"

#include <algorithm>
#include <atomic>
#include <cerrno>
#include <chrono>
#include <cmath>
#include <cstring>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <new>
#include <stdexcept>
#include <string>
#include <vector>

static std::atomic<bool> g_track_allocations(false);
static std::atomic<uint64_t> g_tracked_allocations(0);

#if defined(_MSC_VER)
#define LEO2_BENCH_NOINLINE __declspec(noinline)
#elif defined(__GNUC__) || defined(__clang__)
#define LEO2_BENCH_NOINLINE __attribute__((noinline))
#else
#define LEO2_BENCH_NOINLINE
#endif

LEO2_BENCH_NOINLINE void* operator new(size_t bytes)
{
    if (g_track_allocations.load(std::memory_order_relaxed))
        g_tracked_allocations.fetch_add(1, std::memory_order_relaxed);
    void* const result = malloc(bytes == 0 ? 1 : bytes);
    if (!result)
        throw std::bad_alloc();
    return result;
}

LEO2_BENCH_NOINLINE void* operator new[](size_t bytes)
{
    return ::operator new(bytes);
}

LEO2_BENCH_NOINLINE void* operator new(
    size_t bytes, const std::nothrow_t&) noexcept
{
    try
    {
        return ::operator new(bytes);
    }
    catch (...)
    {
        return NULL;
    }
}

LEO2_BENCH_NOINLINE void* operator new[](
    size_t bytes, const std::nothrow_t&) noexcept
{
    return ::operator new(bytes, std::nothrow);
}

LEO2_BENCH_NOINLINE void operator delete(void* pointer) noexcept
{
    free(pointer);
}

LEO2_BENCH_NOINLINE void operator delete[](void* pointer) noexcept
{
    free(pointer);
}

LEO2_BENCH_NOINLINE void operator delete(
    void* pointer, const std::nothrow_t&) noexcept
{
    ::operator delete(pointer);
}

LEO2_BENCH_NOINLINE void operator delete[](
    void* pointer, const std::nothrow_t&) noexcept
{
    ::operator delete[](pointer);
}

#undef LEO2_BENCH_NOINLINE

namespace {

struct Options
{
    uint32_t k;
    uint32_t r;
    leo2_profile profile;
    leo2_field field;
    leo2_backend backend;
    uint64_t bytes;
    std::string parity;
    size_t iterations;
    size_t warmup;
    size_t codec_repetitions;
    size_t setup_repetitions;
    size_t execute_repetitions;
    std::string output;

    Options()
        : k(240)
        , r(16)
        , profile(LEO2_PROFILE_LEGACY_HIGH_V1)
        , field(LEO2_FIELD_GF8)
        , backend(LEO2_BACKEND_AVX2)
        , bytes(64)
        , parity("missing")
        , iterations(21)
        , warmup(3)
        , codec_repetitions(5)
        , setup_repetitions(10000)
        , execute_repetitions(100000)
        , output("-")
    {}
};

struct Summary
{
    double median_us;
    double mad_us;
    double minimum_us;
    double maximum_us;
    std::vector<double> samples_us;
};

void fail(const std::string& message)
{
    throw std::runtime_error(message);
}

void require_result(leo2_result result, const char* operation)
{
    if (result != LEO2_SUCCESS)
        fail(std::string(operation) + ": " + leo2_result_string(result));
}

uint64_t parse_u64(const char* text, const char* option)
{
    if (!text || !*text || *text == '-')
        fail(std::string("invalid ") + option);
    errno = 0;
    char* end = NULL;
    const unsigned long long value = strtoull(text, &end, 10);
    if (errno == ERANGE || !end || *end != '\0' ||
        value > std::numeric_limits<uint64_t>::max())
        fail(std::string("invalid ") + option);
    return static_cast<uint64_t>(value);
}

uint32_t parse_u32(const char* text, const char* option)
{
    const uint64_t value = parse_u64(text, option);
    if (value > std::numeric_limits<uint32_t>::max())
        fail(std::string("out-of-range ") + option);
    return static_cast<uint32_t>(value);
}

size_t parse_size(const char* text, const char* option)
{
    const uint64_t value = parse_u64(text, option);
    if (value > std::numeric_limits<size_t>::max())
        fail(std::string("out-of-range ") + option);
    return static_cast<size_t>(value);
}

const char* need_value(int argc, char** argv, int& index)
{
    if (++index >= argc)
        fail(std::string("missing value after ") + argv[index - 1]);
    return argv[index];
}

leo2_profile parse_profile(const std::string& text)
{
    if (text == "auto") return LEO2_PROFILE_AUTO;
    if (text == "high") return LEO2_PROFILE_LEGACY_HIGH_V1;
    if (text == "low") return LEO2_PROFILE_LOW_V1;
    fail("invalid --profile");
    return LEO2_PROFILE_AUTO;
}

leo2_field parse_field(const std::string& text)
{
    if (text == "auto") return LEO2_FIELD_AUTO;
    if (text == "gf8") return LEO2_FIELD_GF8;
    if (text == "gf16") return LEO2_FIELD_GF16;
    fail("invalid --field");
    return LEO2_FIELD_AUTO;
}

leo2_backend parse_backend(const std::string& text)
{
    if (text == "auto") return LEO2_BACKEND_AUTO;
    if (text == "scalar") return LEO2_BACKEND_SCALAR;
    if (text == "ssse3") return LEO2_BACKEND_SSSE3;
    if (text == "avx2") return LEO2_BACKEND_AVX2;
    if (text == "avx512") return LEO2_BACKEND_AVX512;
    if (text == "neon") return LEO2_BACKEND_NEON;
    fail("invalid --backend");
    return LEO2_BACKEND_AUTO;
}

Options parse_options(int argc, char** argv)
{
    Options options;
    for (int i = 1; i < argc; ++i)
    {
        const std::string argument(argv[i]);
        if (argument == "--k")
            options.k = parse_u32(need_value(argc, argv, i), "--k");
        else if (argument == "--r")
            options.r = parse_u32(need_value(argc, argv, i), "--r");
        else if (argument == "--profile")
            options.profile = parse_profile(need_value(argc, argv, i));
        else if (argument == "--field")
            options.field = parse_field(need_value(argc, argv, i));
        else if (argument == "--backend")
            options.backend = parse_backend(need_value(argc, argv, i));
        else if (argument == "--bytes")
            options.bytes = parse_u64(need_value(argc, argv, i), "--bytes");
        else if (argument == "--parity")
            options.parity = need_value(argc, argv, i);
        else if (argument == "--iterations")
            options.iterations = parse_size(
                need_value(argc, argv, i), "--iterations");
        else if (argument == "--warmup")
            options.warmup = parse_size(
                need_value(argc, argv, i), "--warmup");
        else if (argument == "--codec-repetitions")
            options.codec_repetitions = parse_size(
                need_value(argc, argv, i), "--codec-repetitions");
        else if (argument == "--setup-repetitions")
            options.setup_repetitions = parse_size(
                need_value(argc, argv, i), "--setup-repetitions");
        else if (argument == "--execute-repetitions")
            options.execute_repetitions = parse_size(
                need_value(argc, argv, i), "--execute-repetitions");
        else if (argument == "--json")
            options.output = need_value(argc, argv, i);
        else if (argument == "--help" || argument == "-h")
        {
            std::cout
                << "Usage: " << argv[0] << " [options]\n"
                << "  --k N --r N --profile auto|high|low\n"
                << "  --field auto|gf8|gf16 --backend auto|scalar|ssse3|avx2|avx512|neon\n"
                << "  --bytes N --parity missing|mixed|all\n"
                << "  --iterations N --warmup N\n"
                << "  --codec-repetitions N --setup-repetitions N\n"
                << "  --execute-repetitions N --json PATH|-\n";
            exit(0);
        }
        else
            fail("unknown option: " + argument);
    }
    if (options.k == 0 || options.r == 0 || options.bytes == 0 ||
        options.iterations == 0 || options.codec_repetitions == 0 ||
        options.setup_repetitions == 0 || options.execute_repetitions == 0)
        fail("counts, bytes, iterations, and repetitions must be positive");
    if (options.parity != "missing" && options.parity != "mixed" &&
        options.parity != "all")
        fail("invalid --parity");
    return options;
}

double median(std::vector<double> values)
{
    std::sort(values.begin(), values.end());
    const size_t middle = values.size() / 2;
    if ((values.size() & 1u) != 0)
        return values[middle];
    return (values[middle - 1] + values[middle]) * 0.5;
}

template <typename Function>
Summary measure(
    size_t iterations,
    size_t warmup,
    size_t repetitions,
    const Function& function)
{
    for (size_t i = 0; i < warmup; ++i)
        function();
    Summary summary;
    summary.samples_us.reserve(iterations);
    for (size_t sample = 0; sample < iterations; ++sample)
    {
        const std::chrono::steady_clock::time_point begin =
            std::chrono::steady_clock::now();
        for (size_t repetition = 0; repetition < repetitions; ++repetition)
            function();
        const std::chrono::steady_clock::time_point end =
            std::chrono::steady_clock::now();
        const double elapsed_us =
            std::chrono::duration_cast<std::chrono::duration<double,
                std::micro> >(end - begin).count();
        summary.samples_us.push_back(
            elapsed_us / static_cast<double>(repetitions));
    }
    summary.median_us = median(summary.samples_us);
    std::vector<double> deviations(summary.samples_us.size());
    for (size_t i = 0; i < summary.samples_us.size(); ++i)
        deviations[i] = std::fabs(summary.samples_us[i] - summary.median_us);
    summary.mad_us = median(deviations);
    summary.minimum_us = *std::min_element(
        summary.samples_us.begin(), summary.samples_us.end());
    summary.maximum_us = *std::max_element(
        summary.samples_us.begin(), summary.samples_us.end());
    return summary;
}

template <typename Function>
uint64_t count_allocations(const Function& function)
{
    struct TrackingGuard
    {
        TrackingGuard()
        {
            g_tracked_allocations.store(0, std::memory_order_relaxed);
            g_track_allocations.store(true, std::memory_order_release);
        }

        ~TrackingGuard()
        {
            g_track_allocations.store(false, std::memory_order_release);
        }
    } guard;
    function();
    return g_tracked_allocations.load(std::memory_order_relaxed);
}

void write_summary(std::ostream& output, const Summary& summary)
{
    output << "{\"median_us\":" << summary.median_us
           << ",\"mad_us\":" << summary.mad_us
           << ",\"minimum_us\":" << summary.minimum_us
           << ",\"maximum_us\":" << summary.maximum_us
           << ",\"samples_us\":[";
    for (size_t i = 0; i < summary.samples_us.size(); ++i)
    {
        if (i != 0)
            output << ',';
        output << summary.samples_us[i];
    }
    output << "]}";
}

const char* profile_name(leo2_profile profile)
{
    if (profile == LEO2_PROFILE_LEGACY_HIGH_V1) return "high";
    if (profile == LEO2_PROFILE_LOW_V1) return "low";
    return "auto";
}

const char* field_name(leo2_field field)
{
    if (field == LEO2_FIELD_GF8) return "gf8";
    if (field == LEO2_FIELD_GF16) return "gf16";
    return "auto";
}

const char* backend_name(leo2_backend backend)
{
    switch (backend)
    {
    case LEO2_BACKEND_SCALAR: return "scalar";
    case LEO2_BACKEND_SSSE3: return "ssse3";
    case LEO2_BACKEND_AVX2: return "avx2";
    case LEO2_BACKEND_AVX512: return "avx512";
    case LEO2_BACKEND_NEON: return "neon";
    default: return "auto";
    }
}

} // namespace

int main(int argc, char** argv)
{
    try
    {
        const Options options = parse_options(argc, argv);
        leo2_context_options context_options;
        memset(&context_options, 0, sizeof(context_options));
        context_options.struct_size = sizeof(context_options);
        context_options.backend = options.backend;
        context_options.thread_count = 1;
        leo2_context* context = NULL;
        require_result(leo2_context_create(&context_options, &context),
            "context create");

        leo2_codec* codec = NULL;
        require_result(leo2_codec_create(context, options.k, options.r,
            options.profile, options.field, NULL, &codec), "codec create");
        std::vector<uint8_t> original_present(options.k, 1);
        std::vector<uint8_t> recovery_present(options.r, 0);
        if (options.parity == "all")
            std::fill(recovery_present.begin(), recovery_present.end(), 1);
        else if (options.parity == "mixed")
            for (uint32_t i = 0; i < options.r; ++i)
                recovery_present[i] = static_cast<uint8_t>(i & 1u);

        leo2_decode_plan* plan = NULL;
        require_result(leo2_decode_plan_create(codec,
            original_present.data(), recovery_present.data(), &plan),
            "persistent plan create");
        require_result(leo2_decode_plan_execute(plan, options.bytes,
            NULL, NULL, NULL, NULL, 0), "persistent no-loss execute");
        require_result(leo2_decode(codec, options.bytes,
            original_present.data(), recovery_present.data(),
            NULL, NULL, NULL, NULL, 0), "one-shot no-loss execute");

        const bool short_circuit = leopard2_internal::
            OneShotNoLossShortCircuitExperimentEnabled();
        const uint64_t plan_setup_allocations = count_allocations([&]() {
            leo2_decode_plan* temporary = NULL;
            require_result(leo2_decode_plan_create(codec,
                original_present.data(), recovery_present.data(), &temporary),
                "allocation-audited plan create");
            leo2_decode_plan_destroy(temporary);
        });
        const uint64_t execute_allocations = count_allocations([&]() {
            require_result(leo2_decode_plan_execute(plan, options.bytes,
                NULL, NULL, NULL, NULL, 0),
                "allocation-audited plan execute");
        });
        const uint64_t plan_first_allocations = count_allocations([&]() {
            leo2_decode_plan* temporary = NULL;
            require_result(leo2_decode_plan_create(codec,
                original_present.data(), recovery_present.data(), &temporary),
                "allocation-audited plan+first create");
            require_result(leo2_decode_plan_execute(temporary, options.bytes,
                NULL, NULL, NULL, NULL, 0),
                "allocation-audited plan+first execute");
            leo2_decode_plan_destroy(temporary);
        });
        const uint64_t one_shot_allocations = count_allocations([&]() {
            require_result(leo2_decode(codec, options.bytes,
                original_present.data(), recovery_present.data(),
                NULL, NULL, NULL, NULL, 0),
                "allocation-audited one-shot execute");
        });
        if (plan_setup_allocations == 0 || plan_first_allocations == 0)
            fail("no-loss plan allocation audit observed no plan allocation");
        if (execute_allocations != 0)
            fail("reusable no-loss plan execution allocated");
        if (short_circuit ? one_shot_allocations != 0
                          : one_shot_allocations == 0)
            fail("one-shot allocation count disagrees with linked marker");

        const Summary codec_setup = measure(options.iterations,
            options.warmup, options.codec_repetitions, [&]() {
                leo2_codec* temporary = NULL;
                require_result(leo2_codec_create(context, options.k,
                    options.r, options.profile, options.field, NULL,
                    &temporary), "timed codec create");
                leo2_codec_destroy(temporary);
            });
        const Summary plan_setup = measure(options.iterations,
            options.warmup, options.setup_repetitions, [&]() {
                leo2_decode_plan* temporary = NULL;
                require_result(leo2_decode_plan_create(codec,
                    original_present.data(), recovery_present.data(),
                    &temporary), "timed plan create");
                leo2_decode_plan_destroy(temporary);
            });
        const Summary plan_execute = measure(options.iterations,
            options.warmup, options.execute_repetitions, [&]() {
                require_result(leo2_decode_plan_execute(plan, options.bytes,
                    NULL, NULL, NULL, NULL, 0), "timed plan execute");
            });
        const Summary plan_first = measure(options.iterations,
            options.warmup, options.setup_repetitions, [&]() {
                leo2_decode_plan* temporary = NULL;
                require_result(leo2_decode_plan_create(codec,
                    original_present.data(), recovery_present.data(),
                    &temporary), "timed plan+first create");
                require_result(leo2_decode_plan_execute(temporary,
                    options.bytes, NULL, NULL, NULL, NULL, 0),
                    "timed plan+first execute");
                leo2_decode_plan_destroy(temporary);
            });
        const Summary one_shot = measure(options.iterations,
            options.warmup, options.setup_repetitions, [&]() {
                require_result(leo2_decode(codec, options.bytes,
                    original_present.data(), recovery_present.data(),
                    NULL, NULL, NULL, NULL, 0), "timed one-shot decode");
            });

        std::ofstream file;
        std::ostream* output = &std::cout;
        if (options.output != "-")
        {
            file.open(options.output.c_str(),
                std::ios::out | std::ios::trunc);
            if (!file)
                fail("failed to open JSON output");
            output = &file;
        }
        *output << std::setprecision(17)
            << "{\n  \"schema\":\"leopard2-no-loss-setup/v1\",\n"
            << "  \"config\":{\"k\":" << options.k
            << ",\"r\":" << options.r
            << ",\"profile_requested\":\""
            << profile_name(options.profile)
            << "\",\"profile_selected\":\""
            << profile_name(leo2_codec_profile(codec))
            << "\",\"field_requested\":\"" << field_name(options.field)
            << "\",\"field_selected\":\""
            << field_name(leo2_codec_field(codec))
            << "\",\"backend_requested\":\""
            << backend_name(options.backend)
            << "\",\"backend_selected\":\""
            << backend_name(leo2_context_backend(context))
            << "\",\"parent_count\":" << leo2_codec_parent_count(codec)
            << ",\"padded_side\":" << leo2_codec_padded_side(codec)
            << ",\"bytes\":" << options.bytes
            << ",\"parity\":\"" << options.parity
            << "\",\"iterations\":" << options.iterations
            << ",\"warmup\":" << options.warmup
            << ",\"codec_repetitions\":" << options.codec_repetitions
            << ",\"setup_repetitions\":" << options.setup_repetitions
            << ",\"execute_repetitions\":"
            << options.execute_repetitions << "},\n"
            << "  \"linked_experiment\":{\"one_shot_no_loss_short_circuit\":"
            << (short_circuit ? "true" : "false")
            << ",\"plan_setup_allocations\":" << plan_setup_allocations
            << ",\"plan_execute_allocations\":" << execute_allocations
            << ",\"plan_plus_first_allocations\":"
            << plan_first_allocations
            << ",\"one_shot_allocations\":" << one_shot_allocations
            << "},\n  \"metrics\":{\n    \"codec_setup\":";
        write_summary(*output, codec_setup);
        *output << ",\n    \"plan_setup\":";
        write_summary(*output, plan_setup);
        *output << ",\n    \"plan_execute\":";
        write_summary(*output, plan_execute);
        *output << ",\n    \"plan_plus_first_execute\":";
        write_summary(*output, plan_first);
        *output << ",\n    \"one_shot_decode\":";
        write_summary(*output, one_shot);
        *output << "\n  },\n  \"reuse\":[";
        static const size_t reuse_counts[] = { 1, 2, 8, 64 };
        for (size_t i = 0;
             i < sizeof(reuse_counts) / sizeof(reuse_counts[0]); ++i)
        {
            if (i != 0)
                *output << ',';
            const double reusable_us =
                (plan_setup.median_us + plan_execute.median_us *
                    static_cast<double>(reuse_counts[i])) /
                static_cast<double>(reuse_counts[i]);
            *output << "{\"count\":" << reuse_counts[i]
                    << ",\"reusable_plan_amortized_us\":" << reusable_us
                    << ",\"one_shot_per_call_us\":"
                    << one_shot.median_us << '}';
        }
        *output << "]\n}\n";
        output->flush();
        if (!*output)
            fail("failed to write JSON output");

        leo2_decode_plan_destroy(plan);
        leo2_codec_destroy(codec);
        leo2_context_destroy(context);
        return 0;
    }
    catch (const std::exception& error)
    {
        g_track_allocations.store(false, std::memory_order_release);
        std::cerr << error.what() << std::endl;
        return 1;
    }
}
