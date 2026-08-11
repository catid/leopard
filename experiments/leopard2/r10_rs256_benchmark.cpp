/*
    Copyright (c) 2017 Christopher A. Taylor.  All rights reserved.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the conditions in License.md
    are met.

    Narrow RS(256,K) decoding benchmark for reproducing the workload in
    Chen et al., "Two Fast Erasure Decoding Algorithms for Reed-Solomon Codes
    Based on LCH-FFT" (R10).  One 1024-byte shard represents the paper's 1024
    codewords that share an erasure-locator set.  Define
    LEO2_R10_LEGACY_ONLY when compiling this adapter against an exact Leopard
    main checkout; otherwise it benchmarks Leopard2's forced Algorithm 4/5 or
    retained generic transform.
*/

#if defined(LEO2_R10_LEGACY_ONLY)
#include "leopard.h"
#else
#include "leopard2.h"
#include "Leopard2Dispatch.h"
#endif

#include <algorithm>
#include <cerrno>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#if defined(_MSC_VER)
#include <malloc.h>
#endif

#ifndef LEO2_R10_SOURCE_COMMIT
#error "LEO2_R10_SOURCE_COMMIT must identify the benchmarked source revision"
#endif

namespace {

static const uint32_t kParent = 256;

enum Mode
{
    ModeAuto,
    ModeSpecialized,
    ModeGeneric,
    ModeLegacy
};

enum Profile
{
    ProfileLow,
    ProfileHigh
};

struct Options
{
    uint32_t k;
    uint64_t bytes;
    uint64_t data_seed;
    uint64_t pattern_seed;
    size_t reuse;
    size_t iterations;
    size_t warmup;
    Mode mode;
    Profile profile;
    std::string output;

    Options()
        : k(240), bytes(1024), data_seed(2), pattern_seed(2), reuse(256),
          iterations(25), warmup(8), mode(
#if defined(LEO2_R10_LEGACY_ONLY)
              ModeLegacy
#else
              ModeSpecialized
#endif
          ), profile(ProfileHigh), output("-")
    {}
};

struct Pattern
{
    std::vector<uint8_t> original_present;
    std::vector<uint8_t> recovery_present;
    std::vector<uint32_t> missing_original;
    std::vector<uint32_t> missing_recovery;
};

struct Summary
{
    double median_us;
    double mad_us;
    double minimum_us;
    double maximum_us;
    std::vector<double> samples_us;
};

class XorShift64
{
public:
    explicit XorShift64(uint64_t seed)
        : state_(seed ? seed : UINT64_C(0x9e3779b97f4a7c15))
    {}

    uint64_t Next()
    {
        uint64_t value = state_;
        value ^= value << 13;
        value ^= value >> 7;
        value ^= value << 17;
        state_ = value;
        return value;
    }

private:
    uint64_t state_;
};

class Fnv1a64
{
public:
    Fnv1a64() : value_(UINT64_C(14695981039346656037)) {}

    void Add(const void* data, size_t bytes)
    {
        const uint8_t* input = static_cast<const uint8_t*>(data);
        for (size_t i = 0; i < bytes; ++i)
        {
            value_ ^= input[i];
            value_ *= UINT64_C(1099511628211);
        }
    }

    std::string Hex() const
    {
        std::ostringstream output;
        output << std::hex << std::setfill('0') << std::setw(16) << value_;
        return output.str();
    }

private:
    uint64_t value_;
};

static void Fail(const std::string& message)
{
    throw std::runtime_error(message);
}

static void* AlignedAllocate(size_t bytes)
{
#if defined(_MSC_VER)
    return _aligned_malloc(bytes, 64);
#else
    void* pointer = NULL;
    return posix_memalign(&pointer, 64, bytes) == 0 ? pointer : NULL;
#endif
}

static void AlignedFree(void* pointer)
{
#if defined(_MSC_VER)
    _aligned_free(pointer);
#else
    free(pointer);
#endif
}

class AlignedBuffer
{
public:
    AlignedBuffer() : data_(NULL), size_(0) {}
    ~AlignedBuffer() { AlignedFree(data_); }

    void Reset(size_t bytes)
    {
        AlignedFree(data_);
        data_ = NULL;
        size_ = 0;
        if (bytes == 0)
            return;
        data_ = AlignedAllocate(bytes);
        if (!data_)
            throw std::bad_alloc();
        size_ = bytes;
        memset(data_, 0, bytes);
    }

    void* data() { return data_; }
    const void* data() const { return data_; }
    uint8_t* bytes() { return static_cast<uint8_t*>(data_); }
    const uint8_t* bytes() const { return static_cast<const uint8_t*>(data_); }
    size_t size() const { return size_; }

private:
    AlignedBuffer(const AlignedBuffer&);
    AlignedBuffer& operator=(const AlignedBuffer&);
    void* data_;
    size_t size_;
};

static uint64_t ParseUnsigned(const std::string& text, const char* option)
{
    if (text.empty() || text[0] == '-')
        Fail(std::string("invalid value for ") + option + ": " + text);
    errno = 0;
    char* end = NULL;
    const unsigned long long value = strtoull(text.c_str(), &end, 10);
    if (errno == ERANGE || end == text.c_str() || *end != '\0')
        Fail(std::string("invalid value for ") + option + ": " + text);
    return static_cast<uint64_t>(value);
}

static size_t ParseSize(const std::string& text, const char* option)
{
    const uint64_t value = ParseUnsigned(text, option);
    if (value > std::numeric_limits<size_t>::max())
        Fail(std::string(option) + " exceeds size_t");
    return static_cast<size_t>(value);
}

static uint32_t ParseUint32(const std::string& text, const char* option)
{
    const uint64_t value = ParseUnsigned(text, option);
    if (value > UINT32_MAX)
        Fail(std::string(option) + " exceeds uint32_t");
    return static_cast<uint32_t>(value);
}

static const char* NeedValue(int argc, char** argv, int& index)
{
    if (++index >= argc)
        Fail(std::string("missing value after ") + argv[index - 1]);
    return argv[index];
}

static void Usage(std::ostream& output, const char* program)
{
    output << "Usage: " << program << " [options]\n"
        << "  --k N                 K in RS(256,K) (default 240)\n"
        << "  --bytes N             Bytes/codewords per erasure group (default 1024)\n"
        << "  --profile low|high    Paper coordinate profile (default from K)\n"
#if !defined(LEO2_R10_LEGACY_ONLY)
        << "  --mode auto|specialized|generic\n"
        << "                         Production AUTO, forced Algorithm 4/5, or generic LCH\n"
#endif
        << "  --pattern-seed N      Deterministic mixed-erasure shuffle seed\n"
        << "  --data-seed N         Deterministic source-data seed\n"
        << "  --reuse N             Calls per sample (default 256)\n"
        << "  --iterations N        Timing samples (default 25)\n"
        << "  --warmup N            Untimed calls (default 8)\n"
        << "  --json PATH           JSON output, or - for stdout\n"
        << "  --help                Show this message\n";
}

static Options ParseOptions(int argc, char** argv)
{
    Options options;
    bool profile_set = false;
    for (int i = 1; i < argc; ++i)
    {
        const std::string argument(argv[i]);
        if (argument == "--help" || argument == "-h")
        {
            Usage(std::cout, argv[0]);
            exit(0);
        }
        else if (argument == "--k")
            options.k = ParseUint32(NeedValue(argc, argv, i), "--k");
        else if (argument == "--bytes")
            options.bytes = ParseUnsigned(NeedValue(argc, argv, i), "--bytes");
        else if (argument == "--pattern-seed")
            options.pattern_seed = ParseUnsigned(
                NeedValue(argc, argv, i), "--pattern-seed");
        else if (argument == "--data-seed")
            options.data_seed = ParseUnsigned(
                NeedValue(argc, argv, i), "--data-seed");
        else if (argument == "--reuse")
            options.reuse = ParseSize(NeedValue(argc, argv, i), "--reuse");
        else if (argument == "--iterations")
            options.iterations = ParseSize(
                NeedValue(argc, argv, i), "--iterations");
        else if (argument == "--warmup")
            options.warmup = ParseSize(NeedValue(argc, argv, i), "--warmup");
        else if (argument == "--profile")
        {
            const std::string profile = NeedValue(argc, argv, i);
            if (profile == "low") options.profile = ProfileLow;
            else if (profile == "high") options.profile = ProfileHigh;
            else Fail("--profile must be low or high");
            profile_set = true;
        }
#if !defined(LEO2_R10_LEGACY_ONLY)
        else if (argument == "--mode")
        {
            const std::string mode = NeedValue(argc, argv, i);
            if (mode == "auto") options.mode = ModeAuto;
            else if (mode == "specialized") options.mode = ModeSpecialized;
            else if (mode == "generic") options.mode = ModeGeneric;
            else Fail("--mode must be auto, specialized, or generic");
        }
#endif
        else if (argument == "--json" || argument == "--output")
            options.output = NeedValue(argc, argv, i);
        else
            Fail("unknown argument: " + argument);
    }

    if (options.k == 0 || options.k >= kParent)
        Fail("--k must be in 1..255");
    if (options.bytes == 0 || options.bytes > std::numeric_limits<size_t>::max())
        Fail("--bytes must be positive and fit size_t");
    if (options.reuse == 0 || options.iterations == 0)
        Fail("--reuse and --iterations must be positive");
    if (!profile_set)
        options.profile = options.k <= 128 ? ProfileLow : ProfileHigh;
#if defined(LEO2_R10_LEGACY_ONLY)
    if (options.profile != ProfileHigh)
        Fail("exact Leopard main only implements the legacy-high profile");
    if (kParent - options.k > options.k)
        Fail("exact Leopard main requires R <= K");
    if ((options.bytes & 63U) != 0)
        Fail("exact Leopard main requires --bytes divisible by 64");
#endif
    return options;
}

static size_t CheckedSize(uint64_t count, uint64_t bytes, const char* what)
{
    if (count != 0 && bytes > std::numeric_limits<size_t>::max() / count)
        Fail(std::string(what) + " size overflow");
    return static_cast<size_t>(count * bytes);
}

static Pattern MakePattern(const Options& options)
{
    const uint32_t r = kParent - options.k;
    std::vector<uint8_t> erased(kParent, 0);
    for (uint32_t i = options.k; i < kParent; ++i)
        erased[i] = 1;
    XorShift64 random(options.pattern_seed ^ UINT64_C(0xd1b54a32d192ed03));
    for (size_t remaining = erased.size(); remaining > 1; --remaining)
    {
        const size_t selected = static_cast<size_t>(random.Next() % remaining);
        std::swap(erased[remaining - 1], erased[selected]);
    }

    Pattern pattern;
    pattern.original_present.assign(options.k, 1);
    pattern.recovery_present.assign(r, 1);
    for (uint32_t coordinate = 0; coordinate < kParent; ++coordinate)
    {
        if (!erased[coordinate])
            continue;
        if (options.profile == ProfileLow)
        {
            if (coordinate < options.k)
            {
                pattern.original_present[coordinate] = 0;
                pattern.missing_original.push_back(coordinate);
            }
            else
            {
                const uint32_t recovery = coordinate - options.k;
                pattern.recovery_present[recovery] = 0;
                pattern.missing_recovery.push_back(recovery);
            }
        }
        else
        {
            if (coordinate < r)
            {
                pattern.recovery_present[coordinate] = 0;
                pattern.missing_recovery.push_back(coordinate);
            }
            else
            {
                const uint32_t original = coordinate - r;
                pattern.original_present[original] = 0;
                pattern.missing_original.push_back(original);
            }
        }
    }
    if (pattern.missing_original.size() + pattern.missing_recovery.size() != r)
        Fail("mixed erasure generator did not produce exactly R erasures");
    if (pattern.missing_original.empty())
        Fail("pattern has no missing originals; choose another --pattern-seed");
    return pattern;
}

static double Median(std::vector<double> values)
{
    const size_t middle = values.size() / 2;
    std::nth_element(values.begin(), values.begin() + middle, values.end());
    const double upper = values[middle];
    if ((values.size() & 1U) != 0)
        return upper;
    std::nth_element(values.begin(), values.begin() + middle - 1,
        values.begin() + middle);
    return (values[middle - 1] + upper) * 0.5;
}

static Summary Summarize(const std::vector<double>& samples)
{
    Summary summary;
    summary.samples_us = samples;
    summary.median_us = Median(samples);
    std::vector<double> deviations(samples.size());
    for (size_t i = 0; i < samples.size(); ++i)
        deviations[i] = std::abs(samples[i] - summary.median_us);
    summary.mad_us = Median(deviations);
    summary.minimum_us = *std::min_element(samples.begin(), samples.end());
    summary.maximum_us = *std::max_element(samples.begin(), samples.end());
    return summary;
}

template <typename Callable>
static Summary Measure(size_t iterations, size_t reuse, Callable callable)
{
    typedef std::chrono::steady_clock Clock;
    std::vector<double> samples;
    samples.reserve(iterations);
    for (size_t sample = 0; sample < iterations; ++sample)
    {
        const Clock::time_point begin = Clock::now();
        for (size_t call = 0; call < reuse; ++call)
            callable();
        const Clock::time_point end = Clock::now();
        const double usec = std::chrono::duration_cast<
            std::chrono::duration<double, std::micro> >(end - begin).count();
        samples.push_back(usec / static_cast<double>(reuse));
    }
    return Summarize(samples);
}

static void WriteIndexArray(
    std::ostream& output, const std::vector<uint32_t>& indices)
{
    output << '[';
    for (size_t i = 0; i < indices.size(); ++i)
    {
        if (i != 0) output << ", ";
        output << indices[i];
    }
    output << ']';
}

static void WriteSummary(
    std::ostream& output, const Summary& summary, uint64_t paper_bytes)
{
    output << "{\n"
        << "      \"median_us\": " << summary.median_us << ",\n"
        << "      \"mad_us\": " << summary.mad_us << ",\n"
        << "      \"minimum_us\": " << summary.minimum_us << ",\n"
        << "      \"maximum_us\": " << summary.maximum_us << ",\n"
        << "      \"paper_throughput_MB_s\": "
        << (static_cast<double>(paper_bytes) / summary.median_us) << ",\n"
        << "      \"samples_us\": [";
    for (size_t i = 0; i < summary.samples_us.size(); ++i)
    {
        if (i != 0) output << ", ";
        output << summary.samples_us[i];
    }
    output << "]\n    }";
}

#if !defined(LEO2_R10_LEGACY_ONLY)

static void RequireLeo2(leo2_result result, const char* operation)
{
    if (result != LEO2_SUCCESS)
    {
        std::ostringstream message;
        message << operation << " failed: " << leo2_result_string(result)
                << " (" << static_cast<int>(result) << ')';
        Fail(message.str());
    }
}

static std::string RunCurrent(const Options& options, const Pattern& pattern)
{
    const uint32_t r = kParent - options.k;
    leo2_context_options context_options = {};
    context_options.struct_size = sizeof(context_options);
    context_options.backend = LEO2_BACKEND_AVX2;
    context_options.thread_count = 1;
    leo2_context* context = NULL;
    RequireLeo2(leo2_context_create(&context_options, &context),
        "context create");

    leo2_codec_options codec_options = {};
    codec_options.struct_size = sizeof(codec_options);
    codec_options.flags = options.mode == ModeGeneric
        ? LEO2_CODEC_FORCE_GENERIC_DECODE
        : options.mode == ModeSpecialized
            ? LEO2_CODEC_FORCE_SPECIALIZED_DECODE |
              LEO2_CODEC_FORCE_TILED_DECODE
            : 0;
    leo2_codec* codec = NULL;
    RequireLeo2(leo2_codec_create(
        context, options.k, r,
        options.profile == ProfileLow
            ? LEO2_PROFILE_LOW_V1 : LEO2_PROFILE_LEGACY_HIGH_V1,
        LEO2_FIELD_GF8, &codec_options, &codec), "codec create");
    if (leo2_context_backend(context) != LEO2_BACKEND_AVX2 ||
        leo2_codec_field(codec) != LEO2_FIELD_GF8 ||
        leo2_codec_parent_count(codec) != kParent)
        Fail("Leopard2 resolved a non-R10 execution identity");

    AlignedBuffer originals;
    AlignedBuffer recovery;
    AlignedBuffer restored;
    originals.Reset(CheckedSize(options.k, options.bytes, "originals"));
    recovery.Reset(CheckedSize(r, options.bytes, "recovery"));
    restored.Reset(CheckedSize(options.k, options.bytes, "restored"));
    XorShift64 data_random(options.data_seed ^ UINT64_C(0x9e3779b97f4a7c15));
    for (size_t i = 0; i < originals.size(); ++i)
        originals.bytes()[i] = static_cast<uint8_t>(data_random.Next() >> 56);

    std::vector<const void*> original(options.k);
    std::vector<const void*> received_original(options.k);
    std::vector<void*> recovery_output(r);
    std::vector<const void*> received_recovery(r);
    std::vector<void*> restored_original(options.k, static_cast<void*>(NULL));
    for (uint32_t i = 0; i < options.k; ++i)
    {
        uint8_t* shard = originals.bytes() +
            static_cast<size_t>(i) * static_cast<size_t>(options.bytes);
        original[i] = shard;
        received_original[i] = pattern.original_present[i] ? shard : NULL;
        if (!pattern.original_present[i])
            restored_original[i] = restored.bytes() +
                static_cast<size_t>(i) * static_cast<size_t>(options.bytes);
    }
    for (uint32_t i = 0; i < r; ++i)
    {
        uint8_t* shard = recovery.bytes() +
            static_cast<size_t>(i) * static_cast<size_t>(options.bytes);
        recovery_output[i] = shard;
        received_recovery[i] = pattern.recovery_present[i] ? shard : NULL;
    }

    size_t encode_scratch_bytes = 0;
    RequireLeo2(leo2_encode_scratch_size(
        codec, options.bytes, &encode_scratch_bytes), "encode scratch query");
    AlignedBuffer encode_scratch;
    encode_scratch.Reset(encode_scratch_bytes);
    RequireLeo2(leo2_encode(codec, options.bytes, &original[0],
        &recovery_output[0], encode_scratch.data(), encode_scratch.size()),
        "encode");

    leo2_decode_plan* plan = NULL;
    RequireLeo2(leo2_decode_plan_create(codec, &pattern.original_present[0],
        &pattern.recovery_present[0], &plan), "plan create");
    if (leo2_decode_plan_missing_original_count(plan) !=
        pattern.missing_original.size())
        Fail("decode plan missing-original count differs from pattern");
    leopard2_internal::DecodePathInfo path_info = {};
    RequireLeo2(leopard2_internal::GetDecodePlanPathInfo(
        plan, options.bytes, false, &path_info), "decode path introspection");
    size_t execution_scratch_bytes = 0;
    size_t one_shot_scratch_bytes = 0;
    RequireLeo2(leo2_decode_plan_scratch_size(
        plan, options.bytes, &execution_scratch_bytes),
        "plan scratch query");
    RequireLeo2(leo2_decode_scratch_size(
        codec, options.bytes, &one_shot_scratch_bytes),
        "one-shot scratch query");
    AlignedBuffer execution_scratch;
    AlignedBuffer one_shot_scratch;
    execution_scratch.Reset(execution_scratch_bytes);
    one_shot_scratch.Reset(one_shot_scratch_bytes);

    const auto run_execution = [&]() {
        RequireLeo2(leo2_decode_plan_execute(plan, options.bytes,
            &received_original[0], &received_recovery[0],
            &restored_original[0], execution_scratch.data(),
            execution_scratch.size()), "plan execute");
    };
    const auto run_one_shot = [&]() {
        RequireLeo2(leo2_decode(codec, options.bytes,
            &pattern.original_present[0], &pattern.recovery_present[0],
            &received_original[0], &received_recovery[0],
            &restored_original[0], one_shot_scratch.data(),
            one_shot_scratch.size()), "one-shot decode");
    };
    const auto check_restored = [&]() {
        for (size_t loss = 0; loss < pattern.missing_original.size(); ++loss)
        {
            const uint32_t index = pattern.missing_original[loss];
            const uint8_t* expected = originals.bytes() +
                static_cast<size_t>(index) * static_cast<size_t>(options.bytes);
            const uint8_t* actual = restored.bytes() +
                static_cast<size_t>(index) * static_cast<size_t>(options.bytes);
            if (memcmp(expected, actual, static_cast<size_t>(options.bytes)) != 0)
                Fail("Leopard2 restored an original incorrectly");
        }
    };

    memset(restored.data(), 0xa5, restored.size());
    run_execution();
    check_restored();
    memset(restored.data(), 0x5a, restored.size());
    run_one_shot();
    check_restored();
    const Summary plan_setup = Measure(options.iterations, 1, [&]() {
        leo2_decode_plan* temporary = NULL;
        RequireLeo2(leo2_decode_plan_create(codec,
            &pattern.original_present[0], &pattern.recovery_present[0],
            &temporary), "timed plan create");
        leo2_decode_plan_destroy(temporary);
    });
    for (size_t i = 0; i < options.warmup; ++i)
        run_execution();
    const Summary execution = Measure(
        options.iterations, options.reuse, run_execution);
    for (size_t i = 0; i < options.warmup; ++i)
        run_one_shot();
    const Summary one_shot = Measure(
        options.iterations, options.reuse, run_one_shot);

    Fnv1a64 original_digest;
    Fnv1a64 parity_digest;
    Fnv1a64 recovered_digest;
    original_digest.Add(originals.data(), originals.size());
    parity_digest.Add(recovery.data(), recovery.size());
    for (size_t loss = 0; loss < pattern.missing_original.size(); ++loss)
    {
        const uint32_t index = pattern.missing_original[loss];
        recovered_digest.Add(restored.bytes() +
            static_cast<size_t>(index) * static_cast<size_t>(options.bytes),
            static_cast<size_t>(options.bytes));
    }

    const uint64_t paper_bytes =
        static_cast<uint64_t>(options.k) * options.bytes;
    std::ostringstream json;
    json << std::fixed << std::setprecision(6)
        << "{\n"
        << "  \"schema\": \"leopard2-r10-rs256-benchmark/v1\",\n"
        << "  \"build\": {\"source_commit\": \""
        << LEO2_R10_SOURCE_COMMIT
        << "\", \"cplusplus\": " << __cplusplus << "},\n"
        << "  \"mode\": \""
        << (options.mode == ModeAuto ? "auto" :
            (options.mode == ModeSpecialized ? "specialized" : "generic"))
        << "\",\n"
        << "  \"algorithm\": \""
        << (options.mode == ModeGeneric ? "generic_lch" :
            (options.profile == ProfileLow ? "algorithm4" : "algorithm5"))
        << "\",\n"
        << "  \"profile\": \""
        << (options.profile == ProfileLow ? "low_v1" : "legacy_high_v1")
        << "\",\n"
        << "  \"field\": \"gf8\",\n"
        << "  \"backend\": \"avx2\",\n"
        << "  \"reusable_decode_path\": \""
        << leopard2_internal::DecodePathName(path_info.path) << "\",\n"
        << "  \"reusable_decode_rule\": \""
        << leopard2_internal::DecodePathRuleName(path_info.rule) << "\",\n"
        << "  \"K\": " << options.k << ",\n"
        << "  \"R\": " << r << ",\n"
        << "  \"N\": 256,\n"
        << "  \"shard_bytes\": " << options.bytes << ",\n"
        << "  \"pattern_seed\": " << options.pattern_seed << ",\n"
        << "  \"data_seed\": " << options.data_seed << ",\n"
        << "  \"missing_original_indices\": ";
    WriteIndexArray(json, pattern.missing_original);
    json << ",\n  \"missing_recovery_indices\": ";
    WriteIndexArray(json, pattern.missing_recovery);
    json << ",\n  \"total_erasure_count\": "
        << pattern.missing_original.size() + pattern.missing_recovery.size()
        << ",\n"
        << "  \"reuse\": " << options.reuse << ",\n"
        << "  \"iterations\": " << options.iterations << ",\n"
        << "  \"correctness\": {\"round_trip\": true},\n"
        << "  \"digests\": {\"algorithm\": \"fnv1a64\", "
        << "\"original_data\": \"" << original_digest.Hex() << "\", "
        << "\"transmitted_parity\": \"" << parity_digest.Hex() << "\", "
        << "\"recovered_originals\": \"" << recovered_digest.Hex()
        << "\"},\n"
        << "  \"scratch\": {\"execution_bytes\": "
        << execution_scratch_bytes << ", \"one_shot_bytes\": "
        << one_shot_scratch_bytes << "},\n"
        << "  \"metrics\": {\n"
        << "    \"plan_setup\": ";
    WriteSummary(json, plan_setup, paper_bytes);
    json << ",\n    \"execution\": ";
    WriteSummary(json, execution, paper_bytes);
    json << ",\n    \"one_shot_including_setup\": ";
    WriteSummary(json, one_shot, paper_bytes);
    json << "\n  }\n}\n";

    leo2_decode_plan_destroy(plan);
    leo2_codec_destroy(codec);
    leo2_context_destroy(context);
    return json.str();
}

#else

static void RequireLegacy(LeopardResult result, const char* operation)
{
    if (result != Leopard_Success)
    {
        std::ostringstream message;
        message << operation << " failed: " << leo_result_string(result)
                << " (" << static_cast<int>(result) << ')';
        Fail(message.str());
    }
}

static std::string RunLegacy(const Options& options, const Pattern& pattern)
{
    const uint32_t r = kParent - options.k;
    RequireLegacy(static_cast<LeopardResult>(leo_init()), "leo_init");
    const unsigned encode_count = leo_encode_work_count(options.k, r);
    const unsigned decode_count = leo_decode_work_count(options.k, r);
    if (encode_count == 0 || decode_count == 0)
        Fail("exact Leopard main rejected the RS(256,K) geometry");

    AlignedBuffer originals;
    AlignedBuffer encode_storage;
    AlignedBuffer decode_storage;
    originals.Reset(CheckedSize(options.k, options.bytes, "originals"));
    encode_storage.Reset(CheckedSize(encode_count, options.bytes, "encode"));
    decode_storage.Reset(CheckedSize(decode_count, options.bytes, "decode"));
    XorShift64 data_random(options.data_seed ^ UINT64_C(0x9e3779b97f4a7c15));
    for (size_t i = 0; i < originals.size(); ++i)
        originals.bytes()[i] = static_cast<uint8_t>(data_random.Next() >> 56);

    std::vector<const void*> original(options.k);
    std::vector<const void*> received_original(options.k);
    std::vector<const void*> received_recovery(r);
    std::vector<void*> encode_work(encode_count);
    std::vector<void*> decode_work(decode_count);
    for (uint32_t i = 0; i < options.k; ++i)
    {
        uint8_t* shard = originals.bytes() +
            static_cast<size_t>(i) * static_cast<size_t>(options.bytes);
        original[i] = shard;
        received_original[i] = pattern.original_present[i] ? shard : NULL;
    }
    for (unsigned i = 0; i < encode_count; ++i)
        encode_work[i] = encode_storage.bytes() +
            static_cast<size_t>(i) * static_cast<size_t>(options.bytes);
    for (unsigned i = 0; i < decode_count; ++i)
        decode_work[i] = decode_storage.bytes() +
            static_cast<size_t>(i) * static_cast<size_t>(options.bytes);
    for (uint32_t i = 0; i < r; ++i)
        received_recovery[i] = pattern.recovery_present[i]
            ? encode_work[i] : NULL;

    RequireLegacy(leo_encode(options.bytes, options.k, r, encode_count,
        &original[0], &encode_work[0]), "encode");
    const auto run_decode = [&]() {
        RequireLegacy(leo_decode(options.bytes, options.k, r, decode_count,
            &received_original[0], &received_recovery[0], &decode_work[0]),
            "decode");
    };
    const auto check_restored = [&]() {
        for (size_t loss = 0; loss < pattern.missing_original.size(); ++loss)
        {
            const uint32_t index = pattern.missing_original[loss];
            const uint8_t* expected = originals.bytes() +
                static_cast<size_t>(index) * static_cast<size_t>(options.bytes);
            const uint8_t* actual = decode_storage.bytes() +
                static_cast<size_t>(index) * static_cast<size_t>(options.bytes);
            if (memcmp(expected, actual, static_cast<size_t>(options.bytes)) != 0)
                Fail("exact Leopard main restored an original incorrectly");
        }
    };
    run_decode();
    check_restored();
    for (size_t i = 0; i < options.warmup; ++i)
        run_decode();
    const Summary decode = Measure(
        options.iterations, options.reuse, run_decode);

    Fnv1a64 original_digest;
    Fnv1a64 parity_digest;
    Fnv1a64 recovered_digest;
    original_digest.Add(originals.data(), originals.size());
    for (uint32_t i = 0; i < r; ++i)
        parity_digest.Add(encode_work[i], static_cast<size_t>(options.bytes));
    for (size_t loss = 0; loss < pattern.missing_original.size(); ++loss)
    {
        const uint32_t index = pattern.missing_original[loss];
        recovered_digest.Add(decode_storage.bytes() +
            static_cast<size_t>(index) * static_cast<size_t>(options.bytes),
            static_cast<size_t>(options.bytes));
    }
    const uint64_t paper_bytes =
        static_cast<uint64_t>(options.k) * options.bytes;
    std::ostringstream json;
    json << std::fixed << std::setprecision(6)
        << "{\n"
        << "  \"schema\": \"leopard2-r10-rs256-benchmark/v1\",\n"
        << "  \"build\": {\"source_commit\": \""
        << LEO2_R10_SOURCE_COMMIT
        << "\", \"cplusplus\": " << __cplusplus << "},\n"
        << "  \"mode\": \"legacy_main\",\n"
        << "  \"algorithm\": \"generic_lch_leopard_main\",\n"
        << "  \"profile\": \"legacy_high_v1\",\n"
        << "  \"field\": \"gf8\",\n"
        << "  \"backend\": \"avx2\",\n"
        << "  \"K\": " << options.k << ",\n"
        << "  \"R\": " << r << ",\n"
        << "  \"N\": 256,\n"
        << "  \"shard_bytes\": " << options.bytes << ",\n"
        << "  \"pattern_seed\": " << options.pattern_seed << ",\n"
        << "  \"data_seed\": " << options.data_seed << ",\n"
        << "  \"missing_original_indices\": ";
    WriteIndexArray(json, pattern.missing_original);
    json << ",\n  \"missing_recovery_indices\": ";
    WriteIndexArray(json, pattern.missing_recovery);
    json << ",\n  \"total_erasure_count\": "
        << pattern.missing_original.size() + pattern.missing_recovery.size()
        << ",\n"
        << "  \"reuse\": " << options.reuse << ",\n"
        << "  \"iterations\": " << options.iterations << ",\n"
        << "  \"correctness\": {\"round_trip\": true},\n"
        << "  \"digests\": {\"algorithm\": \"fnv1a64\", "
        << "\"original_data\": \"" << original_digest.Hex() << "\", "
        << "\"transmitted_parity\": \"" << parity_digest.Hex() << "\", "
        << "\"recovered_originals\": \"" << recovered_digest.Hex()
        << "\"},\n"
        << "  \"metrics\": {\n"
        << "    \"one_shot_including_setup\": ";
    WriteSummary(json, decode, paper_bytes);
    json << "\n  }\n}\n";
    return json.str();
}

#endif

static void WriteOutput(const std::string& path, const std::string& document)
{
    if (path == "-")
    {
        std::cout << document;
        return;
    }
    std::ofstream output(path.c_str(), std::ios::binary | std::ios::trunc);
    if (!output)
        Fail("cannot open output: " + path);
    output << document;
    if (!output)
        Fail("cannot write output: " + path);
}

} // namespace

int main(int argc, char** argv)
{
    try
    {
        const Options options = ParseOptions(argc, argv);
        const Pattern pattern = MakePattern(options);
#if defined(LEO2_R10_LEGACY_ONLY)
        WriteOutput(options.output, RunLegacy(options, pattern));
#else
        WriteOutput(options.output, RunCurrent(options, pattern));
#endif
        return 0;
    }
    catch (const std::exception& error)
    {
        std::cerr << "r10 RS(256,K) benchmark: " << error.what() << std::endl;
        return 1;
    }
}
