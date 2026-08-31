/*
    Copyright (c) 2017 Christopher A. Taylor.  All rights reserved.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the conditions in the Leopard-RS
    BSD license are met.  See the repository LICENSE file.
*/

#include "Leopard2Direct.h"
#include "direct_oracle.h"
#include "leopard2.h"

#if defined(LEO2_ENABLE_TEST_HOOKS)
#error "the sparse-high AUTO benchmark must link the ordinary archive"
#endif
#if !defined(LEO2_HIGH_SPARSE_AUTO_LIBRARY_TEST_HOOKS) || \
    LEO2_HIGH_SPARSE_AUTO_LIBRARY_TEST_HOOKS != 0
#error "the sparse-high AUTO benchmark requires an explicit no-hook marker"
#endif
#if !defined(LEO2_BENCHMARK_SOURCE_ATTESTATION) || \
    !defined(LEO2_BENCHMARK_SOURCE_ATTESTATION_HEADER)
#error "the sparse-high AUTO benchmark requires generated source attestation"
#endif
#include LEO2_BENCHMARK_SOURCE_ATTESTATION_HEADER
#if !defined(LEO2_BENCHMARK_SOURCE_COMMIT) || \
    !defined(LEO2_BENCHMARK_SOURCE_TREE) || \
    !defined(LEO2_BENCHMARK_SOURCE_TRACKED_DIRTY) || \
    !defined(LEO2_BENCHMARK_BUILD_VARIANT) || \
    !defined(LEO2_BENCHMARK_BUILD_TYPE) || \
    !defined(LEO2_BENCHMARK_BUILD_CONFIGURATION_SCHEMA) || \
    !defined(LEO2_BENCHMARK_BUILD_CONFIGURATION_SHA256)
#error "the sparse-high AUTO benchmark requires complete build attestation"
#endif

#include <algorithm>
#include <cerrno>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <locale>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#if defined(__linux__)
#include <sched.h>
#include <unistd.h>
#endif

#if defined(_MSC_VER)
#include <malloc.h>
#endif

namespace {

typedef std::chrono::steady_clock Clock;
typedef std::vector<uint8_t> Bytes;
typedef std::vector<Bytes> Shards;

const uint8_t kCanary = 0xa5;

enum ApiMode
{
    ApiUnspecified,
    ApiOneShot,
    ApiBatch,
    ApiBinding
};

enum PolicyMode
{
    PolicyUnspecified,
    PolicyTablesOffAutoOff,
    PolicyTablesOnAutoOff,
    PolicyTablesOnAutoOn
};

struct Options
{
    uint32_t k;
    uint32_t r;
    uint64_t shard_bytes;
    uint32_t parity_index;
    ApiMode api;
    PolicyMode policy;
    leo2_backend requested_backend;
    bool backend_set;
    size_t batch;
    uint32_t threads;
    size_t iterations;
    size_t setup_iterations;
    size_t calls_per_sample;
    size_t warmups;
    size_t reuse;
    uint64_t seed;
    uint64_t memory_mib;
    std::string output;

    Options()
        : k(2)
        , r(16)
        , shard_bytes(1024)
        , parity_index(0)
        , api(ApiUnspecified)
        , policy(PolicyUnspecified)
        , requested_backend(LEO2_BACKEND_AUTO)
        , backend_set(false)
        , batch(1)
        , threads(1)
        , iterations(5)
        , setup_iterations(5)
        , calls_per_sample(4)
        , warmups(2)
        , reuse(64)
        , seed(UINT64_C(0x48534155544f5631))
        , memory_mib(512)
        , output("-")
    {}
};

struct Summary
{
    std::vector<double> samples_us;
    double median_us;
    double mad_us;
    double minimum_us;
    double maximum_us;
};

void Fail(const std::string& message)
{
    throw std::runtime_error(message);
}

void Require(bool condition, const char* message)
{
    if (!condition)
        Fail(message);
}

void RequireLeo2(leo2_result result, const char* operation)
{
    if (result == LEO2_SUCCESS)
        return;
    std::ostringstream message;
    message << operation << " failed: " << leo2_result_string(result)
            << " (" << static_cast<int>(result) << ')';
    Fail(message.str());
}

uint64_t ParseUnsigned(const std::string& text, const char* option)
{
    if (text.empty() || text[0] == '-')
        Fail(std::string("invalid value for ") + option + ": " + text);
    errno = 0;
    char* end = NULL;
    const unsigned long long value = std::strtoull(text.c_str(), &end, 10);
    if (errno == ERANGE || end == text.c_str() || *end != '\0')
        Fail(std::string("invalid value for ") + option + ": " + text);
    return static_cast<uint64_t>(value);
}

uint32_t ParseUint32(const std::string& text, const char* option)
{
    const uint64_t value = ParseUnsigned(text, option);
    if (value > std::numeric_limits<uint32_t>::max())
        Fail(std::string(option) + " exceeds uint32");
    return static_cast<uint32_t>(value);
}

size_t ParseSize(const std::string& text, const char* option)
{
    const uint64_t value = ParseUnsigned(text, option);
    if (value > std::numeric_limits<size_t>::max())
        Fail(std::string(option) + " exceeds size_t");
    return static_cast<size_t>(value);
}

uint64_t ParseBytes(const std::string& text)
{
    size_t digits = 0;
    while (digits < text.size() && text[digits] >= '0' && text[digits] <= '9')
        ++digits;
    if (digits == 0)
        Fail("invalid --bytes: " + text);
    const uint64_t base = ParseUnsigned(text.substr(0, digits), "--bytes");
    std::string suffix = text.substr(digits);
    for (size_t i = 0; i < suffix.size(); ++i)
        if (suffix[i] >= 'A' && suffix[i] <= 'Z')
            suffix[i] = static_cast<char>(suffix[i] - 'A' + 'a');
    uint64_t multiplier = 1;
    if (suffix.empty() || suffix == "b")
        multiplier = 1;
    else if (suffix == "k" || suffix == "kb")
        multiplier = UINT64_C(1000);
    else if (suffix == "m" || suffix == "mb")
        multiplier = UINT64_C(1000000);
    else if (suffix == "ki" || suffix == "kib")
        multiplier = UINT64_C(1024);
    else if (suffix == "mi" || suffix == "mib")
        multiplier = UINT64_C(1024) * UINT64_C(1024);
    else
        Fail("invalid suffix for --bytes: " + suffix);
    if (base > std::numeric_limits<uint64_t>::max() / multiplier)
        Fail("--bytes overflows uint64");
    return base * multiplier;
}

const char* NeedValue(int argc, char** argv, int& index)
{
    if (++index >= argc)
        Fail(std::string("missing value after ") + argv[index - 1]);
    return argv[index];
}

ApiMode ParseApi(const std::string& text)
{
    if (text == "one-shot" || text == "one_shot") return ApiOneShot;
    if (text == "batch") return ApiBatch;
    if (text == "binding") return ApiBinding;
    Fail("--api must be one-shot, batch, or binding");
    return ApiUnspecified;
}

PolicyMode ParsePolicy(const std::string& text)
{
    if (text == "tables-off-auto-off") return PolicyTablesOffAutoOff;
    if (text == "tables-on-auto-off") return PolicyTablesOnAutoOff;
    if (text == "tables-on-auto-on") return PolicyTablesOnAutoOn;
    Fail("--policy must be tables-off-auto-off, tables-on-auto-off, or tables-on-auto-on");
    return PolicyUnspecified;
}

leo2_backend ParseBackend(const std::string& text)
{
    if (text == "auto") return LEO2_BACKEND_AUTO;
    if (text == "avx2") return LEO2_BACKEND_AVX2;
    Fail("--backend must be auto or avx2");
    return LEO2_BACKEND_AUTO;
}

bool IsSparseHighCampaignTuple(
    uint32_t original_count,
    uint32_t recovery_count,
    uint64_t shard_bytes)
{
    const bool qualified_shape =
        (original_count == 2 || original_count == 3 ||
         original_count == 4 || original_count == 8 ||
         original_count == 12 || original_count == 16) &&
        (recovery_count == 2 || recovery_count == 4 ||
         recovery_count == 8 || recovery_count == 16);
    if (!qualified_shape)
        return false;
    if (shard_bytes == 4096)
        return true;
    const bool boundary_shape =
        (original_count == 2 && recovery_count == 16) ||
        (original_count == 16 && recovery_count == 2);
    return boundary_shape &&
        (shard_bytes == 1024 || shard_bytes == 1088 ||
         shard_bytes == 2048 || shard_bytes == 4032 ||
         shard_bytes == 4160 || shard_bytes == 65536);
}

bool IsSparseHighDirectCandidateTuple(
    uint32_t original_count,
    uint32_t recovery_count,
    uint64_t shard_bytes)
{
    if (recovery_count != 4 && recovery_count != 8 &&
        recovery_count != 16)
        return false;
    if (shard_bytes == 4096)
        return original_count == 2 || original_count == 3 ||
            original_count == 4 || original_count == 8 ||
            original_count == 12 || original_count == 16;
    return original_count == 2 && recovery_count == 16 &&
        (shard_bytes == 1024 || shard_bytes == 1088 ||
         shard_bytes == 2048 || shard_bytes == 4032 ||
         shard_bytes == 4160 || shard_bytes == 65536);
}

const char* ApiName(ApiMode api)
{
    switch (api)
    {
    case ApiOneShot: return "one_shot";
    case ApiBatch: return "batch";
    case ApiBinding: return "binding";
    default: return "unspecified";
    }
}

const char* PolicyName(PolicyMode policy)
{
    switch (policy)
    {
    case PolicyTablesOffAutoOff: return "tables_off_auto_off";
    case PolicyTablesOnAutoOff: return "tables_on_auto_off";
    case PolicyTablesOnAutoOn: return "tables_on_auto_on";
    default: return "unspecified";
    }
}

const char* BackendName(leo2_backend backend)
{
    switch (backend)
    {
    case LEO2_BACKEND_AUTO: return "auto";
    case LEO2_BACKEND_SCALAR: return "scalar";
    case LEO2_BACKEND_SSSE3: return "ssse3";
    case LEO2_BACKEND_AVX2: return "avx2";
    case LEO2_BACKEND_NEON: return "neon";
    case LEO2_BACKEND_AVX512: return "avx512";
    case LEO2_BACKEND_GFNI: return "avx2-gfni";
    }
    return "unknown";
}

const char* CompilerName()
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

const char* CompilerVersion()
{
#if defined(__clang__)
    return __clang_version__;
#elif defined(__GNUC__)
    return __VERSION__;
#elif defined(_MSC_FULL_VER)
#define LEO2_HSA_STRINGIZE_DETAIL(value) #value
#define LEO2_HSA_STRINGIZE(value) LEO2_HSA_STRINGIZE_DETAIL(value)
    return LEO2_HSA_STRINGIZE(_MSC_FULL_VER);
#else
    return "unknown";
#endif
}

void Usage(std::ostream& output, const char* executable)
{
    output
        << "Usage: " << executable << " --api MODE --policy STATE --backend BACKEND [options]\n"
        << "  --k N --r N --bytes N       Cell identity (defaults 2,16,1KiB)\n"
        << "  --parity-index N            Single requested parity row\n"
        << "  --api one-shot|batch|binding\n"
        << "  --batch 1|4|16              One-shot requires one\n"
        << "  --backend auto|avx2         Caller-requested backend identity\n"
        << "  --threads 1|4               Qualified or fenced thread lane\n"
        << "  --policy tables-off-auto-off|tables-on-auto-off|tables-on-auto-on\n"
        << "  --iterations N --setup-iterations N --calls-per-sample N\n"
        << "  --warmups N --reuse 1|8|64 --seed N --memory-mib N\n"
        << "  --json PATH                 '-' writes stdout\n";
}

Options ParseOptions(int argc, char** argv)
{
    Options options;
    for (int i = 1; i < argc; ++i)
    {
        const std::string argument(argv[i]);
        if (argument == "--help" || argument == "-h")
        {
            Usage(std::cout, argv[0]);
            std::exit(0);
        }
        if (argument == "--force-direct" || argument == "--force-transform" ||
            argument == "--auto" || argument == "--mode")
            Fail("forced/test encode modes are forbidden");
        if (argument == "--k")
            options.k = ParseUint32(NeedValue(argc, argv, i), "--k");
        else if (argument == "--r")
            options.r = ParseUint32(NeedValue(argc, argv, i), "--r");
        else if (argument == "--bytes")
            options.shard_bytes = ParseBytes(NeedValue(argc, argv, i));
        else if (argument == "--parity-index")
            options.parity_index = ParseUint32(
                NeedValue(argc, argv, i), "--parity-index");
        else if (argument == "--api")
            options.api = ParseApi(NeedValue(argc, argv, i));
        else if (argument == "--policy")
            options.policy = ParsePolicy(NeedValue(argc, argv, i));
        else if (argument == "--backend")
        {
            options.requested_backend = ParseBackend(NeedValue(argc, argv, i));
            options.backend_set = true;
        }
        else if (argument == "--batch")
            options.batch = ParseSize(NeedValue(argc, argv, i), "--batch");
        else if (argument == "--threads" || argument == "--thread-count")
            options.threads = ParseUint32(NeedValue(argc, argv, i), "--threads");
        else if (argument == "--iterations")
            options.iterations = ParseSize(
                NeedValue(argc, argv, i), "--iterations");
        else if (argument == "--setup-iterations")
            options.setup_iterations = ParseSize(
                NeedValue(argc, argv, i), "--setup-iterations");
        else if (argument == "--calls-per-sample")
            options.calls_per_sample = ParseSize(
                NeedValue(argc, argv, i), "--calls-per-sample");
        else if (argument == "--warmup" || argument == "--warmups")
            options.warmups = ParseSize(NeedValue(argc, argv, i), "--warmups");
        else if (argument == "--reuse")
            options.reuse = ParseSize(NeedValue(argc, argv, i), "--reuse");
        else if (argument == "--seed")
            options.seed = ParseUnsigned(NeedValue(argc, argv, i), "--seed");
        else if (argument == "--memory-mib")
            options.memory_mib = ParseUnsigned(
                NeedValue(argc, argv, i), "--memory-mib");
        else if (argument == "--json" || argument == "--output")
            options.output = NeedValue(argc, argv, i);
        else
            Fail("unknown argument: " + argument);
    }

    if (options.api == ApiUnspecified)
        Fail("--api is required");
    if (options.policy == PolicyUnspecified)
        Fail("--policy is required");
    if (!options.backend_set)
        Fail("--backend is required");
    if (options.k == 0 || options.r == 0 || options.shard_bytes == 0)
        Fail("K, R, and shard bytes must be positive");
    if (!IsSparseHighCampaignTuple(
            options.k, options.r, options.shard_bytes))
        Fail("K, R, and shard bytes are outside the 36-cell sparse-high "
            "campaign envelope (24 candidates plus 12 side-two controls)");
    if (options.shard_bytes > std::numeric_limits<size_t>::max())
        Fail("--bytes must fit in size_t");
    if (options.parity_index >= options.r)
        Fail("--parity-index must be less than R");
    if (options.batch != 1 && options.batch != 4 && options.batch != 16)
        Fail("--batch must be 1, 4, or 16");
    if (options.api == ApiOneShot && options.batch != 1)
        Fail("one-shot API requires --batch 1");
    if (options.threads != 1 && options.threads != 4)
        Fail("--threads must be 1 or 4");
    if (options.reuse != 1 && options.reuse != 8 && options.reuse != 64)
        Fail("--reuse must be 1, 8, or 64");
    if (options.iterations == 0 || options.setup_iterations == 0 ||
        options.calls_per_sample == 0)
        Fail("iteration and calls-per-sample counts must be positive");
    if (options.iterations > 10000 || options.setup_iterations > 10000 ||
        options.calls_per_sample > 100000 || options.warmups > 10000)
        Fail("timing count exceeds the bounded benchmark contract");
    if (options.memory_mib == 0 || options.memory_mib > 65536)
        Fail("--memory-mib must be in [1,65536]");
    return options;
}

uint64_t CheckedProduct(uint64_t a, uint64_t b, const char* name)
{
    if (a != 0 && b > std::numeric_limits<uint64_t>::max() / a)
        Fail(std::string(name) + " overflows uint64");
    return a * b;
}

uint64_t CheckedAdd(uint64_t a, uint64_t b, const char* name)
{
    if (b > std::numeric_limits<uint64_t>::max() - a)
        Fail(std::string(name) + " overflows uint64");
    return a + b;
}

class AlignedBuffer
{
public:
    explicit AlignedBuffer(size_t bytes)
        : value_(NULL)
        , bytes_(bytes)
    {
        if (bytes == 0)
            return;
#if defined(_MSC_VER)
        value_ = _aligned_malloc(bytes, leo2_scratch_alignment());
#else
        if (posix_memalign(&value_, leo2_scratch_alignment(), bytes) != 0)
            value_ = NULL;
#endif
        if (!value_)
            throw std::bad_alloc();
        std::memset(value_, 0, bytes);
    }

    ~AlignedBuffer()
    {
#if defined(_MSC_VER)
        _aligned_free(value_);
#else
        std::free(value_);
#endif
    }

    void* data() const { return value_; }
    size_t size() const { return bytes_; }

private:
    AlignedBuffer(const AlignedBuffer&);
    AlignedBuffer& operator=(const AlignedBuffer&);
    void* value_;
    size_t bytes_;
};

struct ContextDestroy
{
    void operator()(leo2_context* value) const { leo2_context_destroy(value); }
};

struct CodecDestroy
{
    void operator()(leo2_codec* value) const { leo2_codec_destroy(value); }
};

struct BindingDestroy
{
    void operator()(leo2_encode_batch_binding* value) const
    {
        leo2_encode_batch_binding_destroy(value);
    }
};

typedef std::unique_ptr<leo2_context, ContextDestroy> ContextPtr;
typedef std::unique_ptr<leo2_codec, CodecDestroy> CodecPtr;
typedef std::unique_ptr<leo2_encode_batch_binding, BindingDestroy> BindingPtr;

bool PolicyTablesEnabled(PolicyMode policy)
{
    return policy != PolicyTablesOffAutoOff;
}

bool PolicyAutoEnabled(PolicyMode policy)
{
    return policy == PolicyTablesOnAutoOn;
}

ContextPtr CreateContext(const Options& options)
{
    leo2_context_options context_options = {};
    context_options.struct_size = sizeof(context_options);
    context_options.backend = options.requested_backend;
    context_options.thread_count = options.threads;
    leo2_context* raw = NULL;
    RequireLeo2(leo2_context_create(&context_options, &raw), "context create");
    ContextPtr context(raw);
    Require(context.get() != NULL, "context create returned null");
    Require(
        leopard2_internal::SetContextHighSparseDirectEncodePolicyForDiagnostics(
            context.get(), PolicyTablesEnabled(options.policy),
            PolicyAutoEnabled(options.policy)),
        "sparse-high policy is unavailable or rejected");
    if (leo2_context_backend(context.get()) != LEO2_BACKEND_AVX2)
        Fail("effective backend must be avx2 for sparse-high telemetry");
    if (leo2_context_thread_count(context.get()) != options.threads)
        Fail("effective thread count differs from the requested lane");
    return context;
}

leo2_codec* CreateCodecRaw(leo2_context* context, const Options& options)
{
    leo2_codec_options codec_options = {};
    codec_options.struct_size = sizeof(codec_options);
    codec_options.flags = 0;
    codec_options.shard_layout = LEO2_SHARD_LAYOUT_NATIVE_V1;
    leo2_codec* codec = NULL;
    RequireLeo2(leo2_codec_create(
        context, options.k, options.r, LEO2_PROFILE_LEGACY_HIGH_V1,
        LEO2_FIELD_GF8, &codec_options, &codec), "codec create");
    Require(codec != NULL, "codec create returned null");
    return codec;
}

Shards MakeOriginal(uint32_t k, size_t bytes, uint64_t seed)
{
    Shards original(k, Bytes(bytes, 0));
    uint64_t state = seed;
    for (uint32_t shard = 0; shard < k; ++shard)
        for (size_t offset = 0; offset < bytes; ++offset)
        {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            original[shard][offset] = static_cast<uint8_t>(
                state + shard * 29U + offset * 131U);
        }
    return original;
}

Bytes OracleParity(
    const leopard2_test::BinaryField& field,
    const leopard2_test::Matrix& generator,
    const Shards& original,
    uint32_t parity_index)
{
    Bytes parity(original[0].size(), 0);
    const std::vector<leopard2_test::Element>& row =
        generator[original.size() + parity_index];
    for (size_t offset = 0; offset < parity.size(); ++offset)
    {
        leopard2_test::Element value = 0;
        for (size_t source = 0; source < original.size(); ++source)
            value = field.add(value, field.multiply(
                row[source], original[source][offset]));
        parity[offset] = static_cast<uint8_t>(value);
    }
    return parity;
}

struct Stripe
{
    Stripe(
        const Options& options,
        size_t scratch_bytes,
        uint64_t seed,
        const leopard2_test::BinaryField& field,
        const leopard2_test::Matrix& generator)
        : original(MakeOriginal(
              options.k, static_cast<size_t>(options.shard_bytes), seed))
        , recovery(options.r,
              Bytes(static_cast<size_t>(options.shard_bytes), kCanary))
        , original_ptrs(options.k, NULL)
        , recovery_ptrs(options.r, NULL)
        , scratch(scratch_bytes)
        , expected(OracleParity(
              field, generator, original, options.parity_index))
    {
        for (uint32_t i = 0; i < options.k; ++i)
            original_ptrs[i] = original[i].data();
        recovery_ptrs[options.parity_index] =
            recovery[options.parity_index].data();
        item.shard_bytes = options.shard_bytes;
        item.original = original_ptrs.data();
        item.recovery = recovery_ptrs.data();
        item.scratch = scratch.data();
        item.scratch_bytes = scratch.size();
    }

    Shards original;
    Shards recovery;
    std::vector<const void*> original_ptrs;
    std::vector<void*> recovery_ptrs;
    AlignedBuffer scratch;
    Bytes expected;
    leo2_encode_batch_item item;
};

uint64_t Fnv1a64(uint64_t hash, const uint8_t* data, size_t bytes)
{
    for (size_t i = 0; i < bytes; ++i)
    {
        hash ^= data[i];
        hash *= UINT64_C(1099511628211);
    }
    return hash;
}

uint64_t InputChecksum(
    const std::vector<std::unique_ptr<Stripe> >& stripes)
{
    uint64_t hash = UINT64_C(14695981039346656037);
    for (size_t stripe = 0; stripe < stripes.size(); ++stripe)
        for (size_t shard = 0; shard < stripes[stripe]->original.size(); ++shard)
            hash = Fnv1a64(hash, stripes[stripe]->original[shard].data(),
                stripes[stripe]->original[shard].size());
    return hash;
}

uint64_t ParityChecksum(
    const std::vector<std::unique_ptr<Stripe> >& stripes,
    uint32_t parity_index)
{
    uint64_t hash = UINT64_C(14695981039346656037);
    for (size_t stripe = 0; stripe < stripes.size(); ++stripe)
        hash = Fnv1a64(hash,
            stripes[stripe]->recovery[parity_index].data(),
            stripes[stripe]->recovery[parity_index].size());
    return hash;
}

void ResetOutputs(std::vector<std::unique_ptr<Stripe> >& stripes)
{
    for (size_t stripe = 0; stripe < stripes.size(); ++stripe)
        for (size_t recovery = 0;
             recovery < stripes[stripe]->recovery.size(); ++recovery)
            std::fill(stripes[stripe]->recovery[recovery].begin(),
                stripes[stripe]->recovery[recovery].end(), kCanary);
}

void VerifyOutputs(
    const std::vector<std::unique_ptr<Stripe> >& stripes,
    const Options& options)
{
    for (size_t stripe = 0; stripe < stripes.size(); ++stripe)
    {
        const Stripe& value = *stripes[stripe];
        if (value.recovery[options.parity_index] != value.expected)
            Fail("production parity differs from the independent oracle");
        for (uint32_t recovery = 0; recovery < options.r; ++recovery)
        {
            if (recovery == options.parity_index)
                continue;
            if (value.recovery[recovery] !=
                Bytes(static_cast<size_t>(options.shard_bytes), kCanary))
                Fail("an unrequested parity output changed");
        }
    }
}

double Median(std::vector<double> values)
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

Summary Summarize(const std::vector<double>& samples)
{
    Require(!samples.empty(), "cannot summarize an empty sample set");
    Summary summary;
    summary.samples_us = samples;
    summary.median_us = Median(samples);
    std::vector<double> deviations(samples.size());
    for (size_t i = 0; i < samples.size(); ++i)
        deviations[i] = std::fabs(samples[i] - summary.median_us);
    summary.mad_us = Median(deviations);
    summary.minimum_us = *std::min_element(samples.begin(), samples.end());
    summary.maximum_us = *std::max_element(samples.begin(), samples.end());
    return summary;
}

Summary MeasureCodecSetup(leo2_context* context, const Options& options)
{
    for (size_t i = 0; i < options.warmups; ++i)
    {
        CodecPtr codec(CreateCodecRaw(context, options));
    }
    std::vector<double> samples;
    samples.reserve(options.setup_iterations);
    for (size_t i = 0; i < options.setup_iterations; ++i)
    {
        const Clock::time_point begin = Clock::now();
        CodecPtr codec(CreateCodecRaw(context, options));
        const Clock::time_point end = Clock::now();
        samples.push_back(std::chrono::duration_cast<
            std::chrono::duration<double, std::micro> >(end - begin).count());
    }
    return Summarize(samples);
}

BindingPtr CreateBinding(
    const leo2_codec* codec,
    const std::vector<leo2_encode_batch_item>& items)
{
    leo2_encode_batch_binding* raw = NULL;
    RequireLeo2(leo2_encode_batch_binding_create(
        codec, items.data(), items.size(), &raw), "binding create");
    Require(raw != NULL, "binding create returned null");
    return BindingPtr(raw);
}

Summary MeasureBindingSetup(
    const leo2_codec* codec,
    const std::vector<leo2_encode_batch_item>& items,
    const Options& options)
{
    for (size_t i = 0; i < options.warmups; ++i)
    {
        BindingPtr binding = CreateBinding(codec, items);
    }
    std::vector<double> samples;
    samples.reserve(options.setup_iterations);
    for (size_t i = 0; i < options.setup_iterations; ++i)
    {
        const Clock::time_point begin = Clock::now();
        BindingPtr binding = CreateBinding(codec, items);
        const Clock::time_point end = Clock::now();
        samples.push_back(std::chrono::duration_cast<
            std::chrono::duration<double, std::micro> >(end - begin).count());
    }
    return Summarize(samples);
}

void ExecutePublicApi(
    const leo2_codec* codec,
    const std::vector<leo2_encode_batch_item>& items,
    const leo2_encode_batch_binding* binding,
    ApiMode api)
{
    if (api == ApiOneShot)
    {
        const leo2_encode_batch_item& item = items[0];
        RequireLeo2(leo2_encode(codec, item.shard_bytes, item.original,
            item.recovery, item.scratch, item.scratch_bytes),
            "one-shot encode");
    }
    else if (api == ApiBatch)
    {
        RequireLeo2(leo2_encode_batch(codec, items.data(), items.size()),
            "ordinary batch encode");
    }
    else
    {
        Require(binding != NULL, "binding API has no prepared binding");
        RequireLeo2(leo2_encode_batch_binding_execute(binding),
            "binding execute");
    }
}

Summary MeasureExecution(
    const leo2_codec* codec,
    const std::vector<leo2_encode_batch_item>& items,
    const leo2_encode_batch_binding* binding,
    const Options& options)
{
    for (size_t i = 0; i < options.warmups; ++i)
        ExecutePublicApi(codec, items, binding, options.api);
    std::vector<double> samples;
    samples.reserve(options.iterations);
    for (size_t sample = 0; sample < options.iterations; ++sample)
    {
        const Clock::time_point begin = Clock::now();
        for (size_t call = 0; call < options.calls_per_sample; ++call)
            ExecutePublicApi(codec, items, binding, options.api);
        const Clock::time_point end = Clock::now();
        const double elapsed_us = std::chrono::duration_cast<
            std::chrono::duration<double, std::micro> >(end - begin).count();
        samples.push_back(
            elapsed_us / static_cast<double>(options.calls_per_sample));
    }
    return Summarize(samples);
}

std::string JsonEscape(const char* text)
{
    std::ostringstream output;
    output.imbue(std::locale::classic());
    const unsigned char* cursor =
        reinterpret_cast<const unsigned char*>(text);
    for (; *cursor; ++cursor)
    {
        switch (*cursor)
        {
        case '\\': output << "\\\\"; break;
        case '"': output << "\\\""; break;
        case '\b': output << "\\b"; break;
        case '\f': output << "\\f"; break;
        case '\n': output << "\\n"; break;
        case '\r': output << "\\r"; break;
        case '\t': output << "\\t"; break;
        default:
            if (*cursor < 0x20)
                output << "\\u" << std::hex << std::setw(4)
                       << std::setfill('0') << static_cast<unsigned>(*cursor)
                       << std::dec << std::setfill(' ');
            else
                output << static_cast<char>(*cursor);
        }
    }
    return output.str();
}

std::vector<unsigned> RuntimeAffinity()
{
    std::vector<unsigned> result;
#if defined(__linux__)
    cpu_set_t set;
    CPU_ZERO(&set);
    if (sched_getaffinity(0, sizeof(set), &set) != 0)
        Fail("cannot read runtime CPU affinity");
    for (unsigned cpu = 0; cpu < CPU_SETSIZE; ++cpu)
        if (CPU_ISSET(cpu, &set))
            result.push_back(cpu);
#endif
    return result;
}

std::string RuntimeExecutablePath()
{
#if defined(__linux__)
    std::vector<char> buffer(4096, 0);
    const ssize_t length = readlink(
        "/proc/self/exe", buffer.data(), buffer.size() - 1);
    if (length <= 0 || static_cast<size_t>(length) >= buffer.size())
        Fail("cannot resolve runtime executable");
    buffer[static_cast<size_t>(length)] = '\0';
    return std::string(buffer.data());
#else
    return std::string();
#endif
}

std::string Hex64(uint64_t value)
{
    std::ostringstream output;
    output << "0x" << std::hex << std::setw(16) << std::setfill('0') << value;
    return output.str();
}

double GigabytesPerSecond(uint64_t bytes, double microseconds)
{
    if (!(microseconds > 0.0))
        return 0.0;
    return static_cast<double>(bytes) / (microseconds * 1000.0);
}

void WriteSummary(std::ostream& output, const Summary& summary, unsigned indent)
{
    const std::string spaces(indent, ' ');
    output << "{\n" << spaces << "  \"samples_us\": [";
    for (size_t i = 0; i < summary.samples_us.size(); ++i)
    {
        if (i != 0)
            output << ", ";
        output << summary.samples_us[i];
    }
    output << "],\n"
           << spaces << "  \"median_us\": " << summary.median_us << ",\n"
           << spaces << "  \"mad_us\": " << summary.mad_us << ",\n"
           << spaces << "  \"minimum_us\": " << summary.minimum_us << ",\n"
           << spaces << "  \"maximum_us\": " << summary.maximum_us << "\n"
           << spaces << '}';
}

int Run(const Options& options)
{
    if (!leopard2_internal::HighSparseDirectEncodeEnabled())
        Fail("sparse-high table support is not compiled into the production archive");

    const std::vector<unsigned> allowed_cpus = RuntimeAffinity();
    const std::string executable_path = RuntimeExecutablePath();
    ContextPtr context = CreateContext(options);
    CodecPtr codec(CreateCodecRaw(context.get(), options));

    if (leo2_codec_profile(codec.get()) != LEO2_PROFILE_LEGACY_HIGH_V1 ||
        leo2_codec_field(codec.get()) != LEO2_FIELD_GF8 ||
        leo2_codec_shard_layout(codec.get()) != LEO2_SHARD_LAYOUT_NATIVE_V1)
        Fail("resolved codec identity differs from the fixed benchmark contract");

    leopard2_internal::CodecEncodePathInfo path = {};
    Require(leopard2_internal::GetCodecEncodePathInfo(
        codec.get(), options.shard_bytes, 1, &path),
        "cannot query production encode path");
    const size_t expected_rows =
        PolicyTablesEnabled(options.policy) &&
        IsSparseHighDirectCandidateTuple(
            options.k, options.r, options.shard_bytes)
            ? static_cast<size_t>(options.r) : 0;
    if (path.direct_generator_rows != expected_rows)
        Fail("prepared direct rows differ from the selected sparse-table policy");
    if (!PolicyAutoEnabled(options.policy) && path.auto_direct_selected)
        Fail("AUTO-off policy unexpectedly selected direct execution");

    size_t scratch_bytes = 0;
    RequireLeo2(leo2_encode_scratch_size(
        codec.get(), options.shard_bytes, &scratch_bytes),
        "encode scratch query");

    uint64_t estimated_per_stripe = CheckedProduct(
        CheckedAdd(options.k, options.r, "shard count"),
        options.shard_bytes, "shard storage");
    estimated_per_stripe = CheckedAdd(
        estimated_per_stripe, scratch_bytes, "stripe storage");
    estimated_per_stripe = CheckedAdd(
        estimated_per_stripe, options.shard_bytes, "oracle storage");
    const uint64_t estimated_memory = CheckedProduct(
        estimated_per_stripe, options.batch, "batch storage");
    const uint64_t memory_limit = CheckedProduct(
        options.memory_mib, UINT64_C(1024) * UINT64_C(1024),
        "memory limit");
    if (estimated_memory > memory_limit)
        Fail("estimated benchmark storage exceeds --memory-mib");

    const leopard2_test::BinaryField field =
        leopard2_test::make_legacy_gf8();
    const leopard2_test::Matrix generator =
        leopard2_test::direct_systematic_generator(
            field, leopard2_test::make_profile_layout(
                leopard2_test::kLegacyHigh, options.k, options.r));

    std::vector<std::unique_ptr<Stripe> > stripes;
    stripes.reserve(options.batch);
    std::vector<leo2_encode_batch_item> items(options.batch);
    for (size_t i = 0; i < options.batch; ++i)
    {
        stripes.push_back(std::unique_ptr<Stripe>(new Stripe(
            options, scratch_bytes, options.seed + i, field, generator)));
        items[i] = stripes.back()->item;
    }
    const uint64_t input_checksum_before = InputChecksum(stripes);

    BindingPtr binding;
    if (options.api == ApiBinding)
        binding = CreateBinding(codec.get(), items);

    ResetOutputs(stripes);
    Require(
        leopard2_internal::ArmContextHighSparseEncodeRouteWitnessForDiagnostics(
            context.get()),
        "cannot arm sparse-high route witness");
    const size_t witness_public_executions =
        options.api == ApiBinding ? 2 : 1;
    for (size_t i = 0; i < witness_public_executions; ++i)
        ExecutePublicApi(codec.get(), items, binding.get(), options.api);
    leopard2_internal::HighSparseEncodeRouteWitness witness = {};
    Require(
        leopard2_internal::ReadAndDisarmContextHighSparseEncodeRouteWitnessForDiagnostics(
            context.get(), &witness),
        "cannot read and disarm sparse-high route witness");
    const uint64_t expected_witness_calls = CheckedProduct(
        witness_public_executions, options.batch, "witness calls");
    if (path.auto_direct_selected)
    {
        if (witness.direct_calls != expected_witness_calls ||
            witness.transform_calls != 0)
            Fail("actual route witness did not prove direct execution");
    }
    else if (witness.direct_calls != 0 ||
             witness.transform_calls != expected_witness_calls)
    {
        Fail("actual route witness did not prove transform execution");
    }
    VerifyOutputs(stripes, options);
    if (InputChecksum(stripes) != input_checksum_before)
        Fail("qualification execution modified input bytes");

    ResetOutputs(stripes);
    const Summary codec_setup = MeasureCodecSetup(context.get(), options);
    Summary binding_setup;
    const bool has_binding_setup = options.api == ApiBinding;
    if (has_binding_setup)
        binding_setup = MeasureBindingSetup(codec.get(), items, options);
    const Summary execution = MeasureExecution(
        codec.get(), items, binding.get(), options);

    VerifyOutputs(stripes, options);
    const uint64_t input_checksum_after = InputChecksum(stripes);
    if (input_checksum_after != input_checksum_before)
        Fail("timed execution modified input bytes");
    const uint64_t parity_checksum =
        ParityChecksum(stripes, options.parity_index);

    const uint64_t logical_input_bytes = CheckedProduct(CheckedProduct(
        options.k, options.shard_bytes, "logical input bytes"),
        options.batch, "batch logical input bytes");
    const uint64_t requested_output_bytes = CheckedProduct(
        options.shard_bytes, options.batch, "requested output bytes");
    double amortized_us = execution.median_us +
        codec_setup.median_us / static_cast<double>(options.reuse);
    if (has_binding_setup)
        amortized_us +=
            binding_setup.median_us / static_cast<double>(options.reuse);

    std::ostringstream json;
    json.imbue(std::locale::classic());
    json << std::fixed << std::setprecision(6);
    json << "{\n"
         << "  \"schema\": \"leopard2-high-sparse-auto-benchmark-v1\",\n"
         << "  \"authoritative\": false,\n"
         << "  \"authority_note\": \"raw telemetry is non-authoritative; authority requires the pinned paired runner\",\n"
         << "  \"build\": {\n"
         << "    \"compiler\": \"" << CompilerName() << "\",\n"
         << "    \"compiler_version\": \""
         << JsonEscape(CompilerVersion()) << "\",\n"
         << "    \"cplusplus\": " << __cplusplus << ",\n"
         << "    \"backend_variant\": \""
         << JsonEscape(LEO2_BENCHMARK_BUILD_VARIANT) << "\",\n"
         << "    \"build_type\": \""
         << JsonEscape(LEO2_BENCHMARK_BUILD_TYPE) << "\",\n"
         << "    \"build_configuration_schema\": \""
         << LEO2_BENCHMARK_BUILD_CONFIGURATION_SCHEMA << "\",\n"
         << "    \"build_configuration_sha256\": \""
         << LEO2_BENCHMARK_BUILD_CONFIGURATION_SHA256 << "\",\n"
         << "    \"source_commit\": \""
         << LEO2_BENCHMARK_SOURCE_COMMIT << "\",\n"
         << "    \"source_tree\": \""
         << LEO2_BENCHMARK_SOURCE_TREE << "\",\n"
         << "    \"source_tracked_dirty\": "
         << (LEO2_BENCHMARK_SOURCE_TRACKED_DIRTY ? "true" : "false") << ",\n"
         << "    \"library_test_hooks\": false,\n"
         << "    \"high_sparse_tables_compiled\": "
         << (leopard2_internal::HighSparseDirectEncodeEnabled()
                ? "true" : "false") << ",\n"
         << "    \"high_sparse_auto_compiled_default\": "
         << (leopard2_internal::HighSparseDirectEncodeAutoEnabled()
                ? "true" : "false") << "\n"
         << "  },\n"
         << "  \"runtime\": {\n"
#if defined(__linux__)
         << "    \"linux_procfs_affinity_attested\": true,\n"
#else
         << "    \"linux_procfs_affinity_attested\": false,\n"
#endif
         << "    \"executable_path\": \""
         << JsonEscape(executable_path.c_str()) << "\",\n"
         << "    \"allowed_cpus\": [";
    for (size_t i = 0; i < allowed_cpus.size(); ++i)
    {
        if (i != 0)
            json << ", ";
        json << allowed_cpus[i];
    }
    json << "]\n"
         << "  },\n"
         << "  \"parameters\": {\n"
         << "    \"requested_backend\": \""
         << BackendName(options.requested_backend) << "\",\n"
         << "    \"policy\": \"" << PolicyName(options.policy) << "\",\n"
         << "    \"api\": \"" << ApiName(options.api) << "\",\n"
         << "    \"K\": " << options.k << ",\n"
         << "    \"R\": " << options.r << ",\n"
         << "    \"Q\": 1,\n"
         << "    \"parity_index\": " << options.parity_index << ",\n"
         << "    \"profile\": \"legacy_high_v1\",\n"
         << "    \"field\": \"gf8\",\n"
         << "    \"shard_layout\": \"native_v1\",\n"
         << "    \"codec_flags\": 0,\n"
         << "    \"shard_bytes\": " << options.shard_bytes << ",\n"
         << "    \"batch\": " << options.batch << ",\n"
         << "    \"reuse\": " << options.reuse << ",\n"
         << "    \"iterations\": " << options.iterations << ",\n"
         << "    \"setup_iterations\": "
         << options.setup_iterations << ",\n"
         << "    \"calls_per_sample\": "
         << options.calls_per_sample << ",\n"
         << "    \"warmups\": " << options.warmups << ",\n"
         << "    \"requested_thread_count\": " << options.threads << ",\n"
         << "    \"seed\": " << options.seed << ",\n"
         << "    \"memory_mib\": " << options.memory_mib << "\n"
         << "  },\n"
         << "  \"resolved\": {\n"
         << "    \"effective_backend\": \""
         << BackendName(leo2_context_backend(context.get())) << "\",\n"
         << "    \"thread_count\": "
         << leo2_context_thread_count(context.get()) << ",\n"
         << "    \"parent_count\": "
         << leo2_codec_parent_count(codec.get()) << ",\n"
         << "    \"padded_side\": "
         << leo2_codec_padded_side(codec.get()) << ",\n"
         << "    \"direct_generator_rows\": "
         << path.direct_generator_rows << ",\n"
         << "    \"auto_direct_selected\": "
         << (path.auto_direct_selected ? "true" : "false") << ",\n"
         << "    \"selected_route\": \""
         << (path.auto_direct_selected ? "direct" : "transform") << "\"\n"
         << "  },\n"
         << "  \"qualification\": {\n"
         << "    \"route_witness_armed\": true,\n"
         << "    \"witness_public_executions\": "
         << witness_public_executions << ",\n"
         << "    \"expected_item_calls\": "
         << expected_witness_calls << ",\n"
         << "    \"direct_calls\": " << witness.direct_calls << ",\n"
         << "    \"transform_calls\": " << witness.transform_calls << ",\n"
         << "    \"witness_disabled_before_timing\": true\n"
         << "  },\n"
         << "  \"correctness\": {\n"
         << "    \"oracle_algorithm\": \"legacy-gf8-direct-systematic-generator-v1\",\n"
         << "    \"independent_oracle_match\": true,\n"
         << "    \"input_immutable\": true,\n"
         << "    \"unrequested_outputs_untouched\": true,\n"
         << "    \"post_timing_recheck_match\": true,\n"
         << "    \"input_checksum_fnv1a64\": \""
         << Hex64(input_checksum_after) << "\",\n"
         << "    \"parity_checksum_fnv1a64\": \""
         << Hex64(parity_checksum) << "\"\n"
         << "  },\n"
         << "  \"memory\": {\n"
         << "    \"scratch_alignment\": " << leo2_scratch_alignment() << ",\n"
         << "    \"scratch_bytes_per_item\": " << scratch_bytes << ",\n"
         << "    \"scratch_bytes_batch\": "
         << CheckedProduct(scratch_bytes, options.batch, "batch scratch")
         << ",\n"
         << "    \"logical_input_bytes_per_call\": "
         << logical_input_bytes << ",\n"
         << "    \"requested_output_bytes_per_call\": "
         << requested_output_bytes << ",\n"
         << "    \"estimated_benchmark_storage_bytes\": "
         << estimated_memory << "\n"
         << "  },\n"
         << "  \"metrics\": {\n"
         << "    \"codec_setup\": ";
    WriteSummary(json, codec_setup, 4);
    json << ",\n"
         << "    \"binding_setup\": ";
    if (has_binding_setup)
        WriteSummary(json, binding_setup, 4);
    else
        json << "null";
    json << ",\n"
         << "    \"execution\": ";
    WriteSummary(json, execution, 4);
    json << ",\n"
         << "    \"amortized\": {\n"
         << "      \"reuse_count\": " << options.reuse << ",\n"
         << "      \"derived_median_us_per_api_call\": "
         << amortized_us << ",\n"
         << "      \"logical_input_GB_per_s\": "
         << GigabytesPerSecond(logical_input_bytes, amortized_us) << ",\n"
         << "      \"requested_parity_output_GB_per_s\": "
         << GigabytesPerSecond(requested_output_bytes, amortized_us) << "\n"
         << "    }\n"
         << "  },\n"
         << "  \"methodology\": {\n"
         << "    \"codec_setup_scope\": \"codec_create only; context and policy reused; destroy excluded\",\n"
         << "    \"binding_setup_scope\": \"binding_create only against preallocated descriptors; destroy excluded; null outside binding API\",\n"
         << "    \"execution_scope\": \"one public API call including ordinary validation; calls averaged within each sample\",\n"
         << "    \"route_witness_scope\": \"one untimed public execution, or two for binding reuse proof; disarmed before warmup\",\n"
         << "    \"timing_allocation_scope\": \"all benchmark-owned shards, descriptors, scratch, bindings, and sample vectors are prepared outside execution spans\",\n"
         << "    \"amortization_formula\": \"execution + codec_setup/reuse + binding_setup/reuse when applicable; one binding per codec\",\n"
         << "    \"affinity_scope\": \"affinity is established externally and captured before context creation\",\n"
         << "    \"production_autotuning\": false\n"
         << "  }\n"
         << "}\n";

    if (options.output == "-")
        std::cout << json.str();
    else
    {
        std::ofstream output(options.output.c_str(),
            std::ios::out | std::ios::trunc);
        if (!output)
            Fail("cannot open JSON output: " + options.output);
        output << json.str();
        if (!output)
            Fail("failed writing JSON output: " + options.output);
    }
    return 0;
}

} // namespace

int main(int argc, char** argv)
{
    try
    {
        return Run(ParseOptions(argc, argv));
    }
    catch (const std::exception& error)
    {
        std::cerr << "leopard2 sparse-high AUTO benchmark: "
                  << error.what() << std::endl;
        return 1;
    }
}
