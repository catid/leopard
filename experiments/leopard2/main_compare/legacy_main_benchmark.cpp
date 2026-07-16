/*
    Copyright (c) 2017 Christopher A. Taylor.  All rights reserved.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the conditions in License.md
    are met.

    This adapter intentionally links only the exact Leopard master baseline.
    It mirrors bench/leopard2/benchmark.cpp's deterministic data, loss
    selection, timing summaries, and throughput definitions so separate
    process-level comparisons can use equivalent work.
*/

#include "leopard.h"

#include <algorithm>
#include <cerrno>
#include <chrono>
#include <cctype>
#include <cmath>
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

#if defined(_MSC_VER)
#include <malloc.h>
#endif

#if defined(_OPENMP)
#include <omp.h>
#endif

#ifndef LEOPARD_MAIN_SOURCE_COMMIT
#error "LEOPARD_MAIN_SOURCE_COMMIT must identify the verified source revision"
#endif

namespace {

struct Options
{
    unsigned k;
    unsigned r;
    uint64_t bytes;
    unsigned losses;
    size_t batch;
    size_t reuse;
    size_t iterations;
    size_t warmup;
    unsigned threads;
    uint64_t seed;
    std::string output;

    Options()
        : k(240), r(16), bytes(65536), losses(1), batch(1), reuse(8),
          iterations(9), warmup(2), threads(1), seed(1), output("-")
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

    uint8_t* bytes() { return static_cast<uint8_t*>(data_); }
    size_t size() const { return size_; }

private:
    AlignedBuffer(const AlignedBuffer&);
    AlignedBuffer& operator=(const AlignedBuffer&);
    void* data_;
    size_t size_;
};

struct Stripe
{
    AlignedBuffer original_storage;
    AlignedBuffer encode_storage;
    AlignedBuffer decode_storage;
    std::vector<const void*> original;
    std::vector<const void*> received_original;
    std::vector<const void*> received_recovery;
    std::vector<void*> encode_work;
    std::vector<void*> decode_work;
};

static size_t CheckedSize(uint64_t count, uint64_t bytes, const char* what)
{
    if (count != 0 && bytes > std::numeric_limits<size_t>::max() / count)
        Fail(std::string(what) + " size overflow");
    return static_cast<size_t>(count * bytes);
}

static uint64_t CheckedProduct(uint64_t left, uint64_t right, const char* what)
{
    if (left != 0 && right > std::numeric_limits<uint64_t>::max() / left)
        Fail(std::string(what) + " overflow");
    return left * right;
}

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

static unsigned ParseUnsignedInt(const std::string& text, const char* option)
{
    const uint64_t value = ParseUnsigned(text, option);
    if (value > std::numeric_limits<unsigned>::max())
        Fail(std::string("value for ") + option + " exceeds unsigned");
    return static_cast<unsigned>(value);
}

static size_t ParseSize(const std::string& text, const char* option)
{
    const uint64_t value = ParseUnsigned(text, option);
    if (value > std::numeric_limits<size_t>::max())
        Fail(std::string("value for ") + option + " exceeds size_t");
    return static_cast<size_t>(value);
}

static uint64_t ParseBytes(const std::string& text)
{
    size_t digits = 0;
    while (digits < text.size() && text[digits] >= '0' && text[digits] <= '9')
        ++digits;
    if (digits == 0)
        Fail("invalid value for --bytes: " + text);
    const uint64_t base = ParseUnsigned(text.substr(0, digits), "--bytes");
    std::string suffix = text.substr(digits);
    std::transform(suffix.begin(), suffix.end(), suffix.begin(),
        [](unsigned char character) {
            return static_cast<char>(std::tolower(character));
        });
    uint64_t multiplier = 1;
    if (suffix.empty() || suffix == "b") multiplier = 1;
    else if (suffix == "k" || suffix == "kb") multiplier = 1000;
    else if (suffix == "m" || suffix == "mb") multiplier = UINT64_C(1000000);
    else if (suffix == "g" || suffix == "gb") multiplier = UINT64_C(1000000000);
    else if (suffix == "ki" || suffix == "kib") multiplier = 1024;
    else if (suffix == "mi" || suffix == "mib") multiplier = UINT64_C(1024) * 1024;
    else if (suffix == "gi" || suffix == "gib") multiplier = UINT64_C(1024) * 1024 * 1024;
    else Fail("invalid suffix for --bytes: " + suffix);
    if (base > std::numeric_limits<uint64_t>::max() / multiplier)
        Fail("--bytes overflow");
    return base * multiplier;
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
        << "  --k N          Original shard count (default 240)\n"
        << "  --r N          Recovery shard count (default 16)\n"
        << "  --bytes N      Bytes per shard; B/kB/MB/GiB suffixes accepted\n"
        << "  --loss N       Missing original shards (default 1)\n"
        << "  --batch N      Independent stripes per timed call (default 1)\n"
        << "  --reuse N      Calls per timing sample (default 8)\n"
        << "  --iterations N Timing samples (default 9)\n"
        << "  --warmup N     Untimed calls (default 2)\n"
        << "  --threads N    OpenMP thread count (default 1)\n"
        << "  --seed N       Deterministic seed (default 1)\n"
        << "  --json PATH    JSON output path, or - for stdout\n"
        << "  --help          Show this message\n";
}

static Options ParseOptions(int argc, char** argv)
{
    Options options;
    for (int i = 1; i < argc; ++i)
    {
        const std::string argument(argv[i]);
        if (argument == "--help" || argument == "-h")
        {
            Usage(std::cout, argv[0]);
            exit(0);
        }
        else if (argument == "--k") options.k = ParseUnsignedInt(NeedValue(argc, argv, i), "--k");
        else if (argument == "--r") options.r = ParseUnsignedInt(NeedValue(argc, argv, i), "--r");
        else if (argument == "--bytes") options.bytes = ParseBytes(NeedValue(argc, argv, i));
        else if (argument == "--loss" || argument == "--losses") options.losses = ParseUnsignedInt(NeedValue(argc, argv, i), "--loss");
        else if (argument == "--batch") options.batch = ParseSize(NeedValue(argc, argv, i), "--batch");
        else if (argument == "--reuse") options.reuse = ParseSize(NeedValue(argc, argv, i), "--reuse");
        else if (argument == "--iterations") options.iterations = ParseSize(NeedValue(argc, argv, i), "--iterations");
        else if (argument == "--warmup") options.warmup = ParseSize(NeedValue(argc, argv, i), "--warmup");
        else if (argument == "--threads") options.threads = ParseUnsignedInt(NeedValue(argc, argv, i), "--threads");
        else if (argument == "--seed") options.seed = ParseUnsigned(NeedValue(argc, argv, i), "--seed");
        else if (argument == "--json" || argument == "--output") options.output = NeedValue(argc, argv, i);
        else Fail("unknown argument: " + argument);
    }
    if (options.k == 0 || options.r == 0 || options.r > options.k)
        Fail("exact Leopard main requires positive counts with R <= K");
    if (options.bytes == 0 || (options.bytes & 63u) != 0 ||
        options.bytes > std::numeric_limits<size_t>::max())
        Fail("exact Leopard main requires positive shard bytes divisible by 64");
    if (options.losses > options.k || options.losses > options.r)
        Fail("--loss cannot exceed K or R");
    if (options.batch == 0 || options.reuse == 0 ||
        options.iterations == 0 || options.threads == 0)
        Fail("--batch, --reuse, --iterations, and --threads must be positive");
    return options;
}

static std::vector<unsigned> SelectLosses(const Options& options)
{
    std::vector<unsigned> order(options.k);
    for (unsigned i = 0; i < options.k; ++i)
        order[i] = i;
    XorShift64 random(options.seed ^ UINT64_C(0xd1b54a32d192ed03));
    for (size_t remaining = order.size(); remaining > 1; --remaining)
    {
        const size_t selected = static_cast<size_t>(random.Next() % remaining);
        std::swap(order[remaining - 1], order[selected]);
    }
    order.resize(options.losses);
    std::sort(order.begin(), order.end());
    return order;
}

static bool IsLost(const std::vector<unsigned>& losses, unsigned index)
{
    return std::binary_search(losses.begin(), losses.end(), index);
}

static void FillOriginals(Stripe& stripe, const Options& options, size_t stripe_index)
{
    XorShift64 random(options.seed ^
        (UINT64_C(0x9e3779b97f4a7c15) * (stripe_index + 1)));
    const size_t total = CheckedSize(options.k, options.bytes, "original data");
    for (size_t i = 0; i < total; ++i)
        stripe.original_storage.bytes()[i] = static_cast<uint8_t>(random.Next() >> 56);
}

static void RequireLegacy(LeopardResult result, const char* operation)
{
    if (result != Leopard_Success)
    {
        std::ostringstream stream;
        stream << operation << " failed: " << leo_result_string(result)
               << " (" << static_cast<int>(result) << ')';
        Fail(stream.str());
    }
}

static void InitializeStripe(Stripe& stripe, const Options& options,
    const std::vector<unsigned>& losses, size_t stripe_index,
    unsigned encode_count, unsigned decode_count)
{
    stripe.original_storage.Reset(CheckedSize(options.k, options.bytes, "original"));
    stripe.encode_storage.Reset(CheckedSize(encode_count, options.bytes, "encode work"));
    stripe.decode_storage.Reset(CheckedSize(decode_count, options.bytes, "decode work"));
    stripe.original.resize(options.k);
    stripe.received_original.resize(options.k);
    stripe.received_recovery.resize(options.r);
    stripe.encode_work.resize(encode_count);
    stripe.decode_work.resize(decode_count);
    FillOriginals(stripe, options, stripe_index);
    for (unsigned i = 0; i < options.k; ++i)
    {
        uint8_t* data = stripe.original_storage.bytes() +
            static_cast<size_t>(i) * static_cast<size_t>(options.bytes);
        stripe.original[i] = data;
        stripe.received_original[i] = IsLost(losses, i) ? NULL : data;
    }
    for (unsigned i = 0; i < encode_count; ++i)
        stripe.encode_work[i] = stripe.encode_storage.bytes() +
            static_cast<size_t>(i) * static_cast<size_t>(options.bytes);
    for (unsigned i = 0; i < decode_count; ++i)
        stripe.decode_work[i] = stripe.decode_storage.bytes() +
            static_cast<size_t>(i) * static_cast<size_t>(options.bytes);
    for (unsigned i = 0; i < options.r; ++i)
        stripe.received_recovery[i] = stripe.encode_work[i];
}

static double Median(std::vector<double> values)
{
    const size_t middle = values.size() / 2;
    std::nth_element(values.begin(), values.begin() + middle, values.end());
    const double upper = values[middle];
    if ((values.size() & 1u) != 0)
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
        deviations[i] = std::fabs(samples[i] - summary.median_us);
    summary.mad_us = Median(deviations);
    summary.minimum_us = *std::min_element(samples.begin(), samples.end());
    summary.maximum_us = *std::max_element(samples.begin(), samples.end());
    return summary;
}

template<class Callable>
static Summary Measure(size_t iterations, size_t inner_calls,
    const Callable& callable)
{
    typedef std::chrono::steady_clock Clock;
    std::vector<double> samples;
    samples.reserve(iterations);
    for (size_t sample = 0; sample < iterations; ++sample)
    {
        const Clock::time_point begin = Clock::now();
        for (size_t call = 0; call < inner_calls; ++call)
            callable();
        const Clock::time_point end = Clock::now();
        const double microseconds = std::chrono::duration_cast<
            std::chrono::duration<double, std::micro> >(end - begin).count();
        samples.push_back(microseconds / static_cast<double>(inner_calls));
    }
    return Summarize(samples);
}

static double Rate(uint64_t bytes, double microseconds)
{
    return bytes == 0 ? 0.0 : static_cast<double>(bytes) / (microseconds * 1000.0);
}

static void WriteRate(std::ostream& output, uint64_t bytes, double microseconds)
{
    if (bytes == 0)
        output << "null";
    else
        output << Rate(bytes, microseconds);
}

static void WriteSummary(std::ostream& output, const Summary& summary,
    uint64_t input_bytes, const char* input_name,
    uint64_t output_bytes, const char* output_name, unsigned indent)
{
    const std::string spaces(indent, ' ');
    output << "{\n"
        << spaces << "  \"median_us_per_batch_call\": " << summary.median_us << ",\n"
        << spaces << "  \"mad_us_per_batch_call\": " << summary.mad_us << ",\n"
        << spaces << "  \"minimum_us_per_batch_call\": " << summary.minimum_us << ",\n"
        << spaces << "  \"maximum_us_per_batch_call\": " << summary.maximum_us << ",\n"
        << spaces << "  \"samples_us_per_batch_call\": [";
    for (size_t i = 0; i < summary.samples_us.size(); ++i)
    {
        if (i != 0) output << ", ";
        output << summary.samples_us[i];
    }
    output << "],\n"
        << spaces << "  \"" << input_name << "\": ";
    WriteRate(output, input_bytes, summary.median_us);
    output << ",\n" << spaces << "  \"" << output_name << "\": ";
    WriteRate(output, output_bytes, summary.median_us);
    output << "\n" << spaces << '}';
}

struct Geometry
{
    unsigned padded_r;
    unsigned parent;
    unsigned encode_work;
    unsigned decode_work;
};

static bool NextPowerOfTwo(unsigned value, unsigned& result)
{
    result = 1;
    while (result < value)
    {
        if (result > std::numeric_limits<unsigned>::max() / 2)
            return false;
        result <<= 1;
    }
    return true;
}

static Geometry ValidateGeometry(const Options& options)
{
    Geometry geometry;
    if (!NextPowerOfTwo(options.r, geometry.padded_r) ||
        geometry.padded_r > std::numeric_limits<unsigned>::max() - options.k ||
        !NextPowerOfTwo(options.k + geometry.padded_r, geometry.parent) ||
        geometry.parent > 65536)
        Fail("requested counts overflow or exceed Leopard main's GF16 parent");
    if (options.k == 1)
        geometry.encode_work = options.r;
    else if (options.r == 1)
        geometry.encode_work = 1;
    else
    {
        if (geometry.padded_r > std::numeric_limits<unsigned>::max() / 2)
            Fail("encode work count overflow");
        geometry.encode_work = geometry.padded_r * 2;
    }
    geometry.decode_work = (options.k == 1 || options.r == 1)
        ? options.k : geometry.parent;
    return geometry;
}

static std::string Run(const Options& options)
{
    unsigned resolved_threads = 1;
#if defined(_OPENMP)
    omp_set_dynamic(0);
    omp_set_num_threads(static_cast<int>(options.threads));
    resolved_threads = static_cast<unsigned>(omp_get_max_threads());
    if (resolved_threads != options.threads)
        Fail("OpenMP did not accept the requested --threads value");
#else
    if (options.threads != 1)
        Fail("this exact main build has no OpenMP support; --threads must be 1");
#endif

    const Geometry geometry = ValidateGeometry(options);
    RequireLegacy(static_cast<LeopardResult>(leo_init()), "leo_init");
    const unsigned encode_count = leo_encode_work_count(options.k, options.r);
    const unsigned decode_count = leo_decode_work_count(options.k, options.r);
    if (encode_count == 0 || decode_count == 0 ||
        encode_count != geometry.encode_work ||
        decode_count != geometry.decode_work)
        Fail("exact Leopard main work counts disagree with validated geometry");

    const std::vector<unsigned> losses = SelectLosses(options);
    std::vector<std::unique_ptr<Stripe> > stripes;
    stripes.reserve(options.batch);
    for (size_t i = 0; i < options.batch; ++i)
    {
        stripes.push_back(std::unique_ptr<Stripe>(new Stripe));
        InitializeStripe(*stripes.back(), options, losses, i,
            encode_count, decode_count);
    }

    // Correctness and workload fingerprints are deliberately outside all
    // timed regions.  Digest traversal order is stripe, shard, then byte.
    Fnv1a64 original_digest;
    Fnv1a64 parity_digest;
    Fnv1a64 recovered_digest;
    for (size_t i = 0; i < stripes.size(); ++i)
    {
        Stripe& stripe = *stripes[i];
        for (unsigned original = 0; original < options.k; ++original)
            original_digest.Add(stripe.original[original],
                static_cast<size_t>(options.bytes));
        RequireLegacy(leo_encode(options.bytes, options.k, options.r,
            encode_count, &stripe.original[0], &stripe.encode_work[0]),
            "correctness encode");
        for (unsigned parity = 0; parity < options.r; ++parity)
            parity_digest.Add(stripe.encode_work[parity],
                static_cast<size_t>(options.bytes));
        RequireLegacy(leo_decode(options.bytes, options.k, options.r,
            decode_count, &stripe.received_original[0],
            &stripe.received_recovery[0], &stripe.decode_work[0]),
            "correctness decode");
        for (size_t loss_i = 0; loss_i < losses.size(); ++loss_i)
        {
            const unsigned index = losses[loss_i];
            const uint8_t* expected = stripe.original_storage.bytes() +
                static_cast<size_t>(index) * static_cast<size_t>(options.bytes);
            const uint8_t* actual = stripe.decode_storage.bytes() +
                static_cast<size_t>(index) * static_cast<size_t>(options.bytes);
            if (memcmp(expected, actual, static_cast<size_t>(options.bytes)) != 0)
                Fail("exact Leopard main restored a shard incorrectly");
            recovered_digest.Add(actual, static_cast<size_t>(options.bytes));
        }
    }

    for (size_t warmup = 0; warmup < options.warmup; ++warmup)
    {
        for (size_t i = 0; i < stripes.size(); ++i)
        {
            Stripe& stripe = *stripes[i];
            RequireLegacy(leo_encode(options.bytes, options.k, options.r,
                encode_count, &stripe.original[0], &stripe.encode_work[0]),
                "encode warmup");
            RequireLegacy(leo_decode(options.bytes, options.k, options.r,
                decode_count, &stripe.received_original[0],
                &stripe.received_recovery[0], &stripe.decode_work[0]),
                "decode warmup");
        }
    }

    const Summary encode = Measure(options.iterations, options.reuse, [&]() {
        for (size_t i = 0; i < stripes.size(); ++i)
        {
            Stripe& stripe = *stripes[i];
            RequireLegacy(leo_encode(options.bytes, options.k, options.r,
                encode_count, &stripe.original[0], &stripe.encode_work[0]),
                "timed encode");
        }
    });
    const Summary decode = Measure(options.iterations, options.reuse, [&]() {
        for (size_t i = 0; i < stripes.size(); ++i)
        {
            Stripe& stripe = *stripes[i];
            RequireLegacy(leo_decode(options.bytes, options.k, options.r,
                decode_count, &stripe.received_original[0],
                &stripe.received_recovery[0], &stripe.decode_work[0]),
                "timed decode");
        }
    });

    const uint64_t batch = static_cast<uint64_t>(options.batch);
    const uint64_t encode_input = CheckedProduct(
        CheckedProduct(options.k, options.bytes, "encode input"), batch,
        "encode input");
    const uint64_t encode_output = CheckedProduct(
        CheckedProduct(options.r, options.bytes, "encode output"), batch,
        "encode output");
    const uint64_t decode_input_shards = options.k - options.losses + options.r;
    const uint64_t decode_input = CheckedProduct(
        CheckedProduct(decode_input_shards, options.bytes, "decode input"),
        batch, "decode input");
    const uint64_t decode_output = CheckedProduct(
        CheckedProduct(options.losses, options.bytes, "decode output"), batch,
        "decode output");
    std::ostringstream json;
    json << std::fixed << std::setprecision(6)
        << "{\n"
        << "  \"schema\": \"leopard-main-benchmark-v1\",\n"
        << "  \"build\": {\n"
        << "    \"main_source_commit\": \"" << LEOPARD_MAIN_SOURCE_COMMIT << "\",\n"
        << "    \"cplusplus\": " << __cplusplus << "\n"
        << "  },\n"
        << "  \"parameters\": {\n"
        << "    \"K\": " << options.k << ",\n"
        << "    \"R\": " << options.r << ",\n"
        << "    \"shard_bytes\": " << options.bytes << ",\n"
        << "    \"loss_count\": " << options.losses << ",\n"
        << "    \"missing_original_indices\": [";
    for (size_t i = 0; i < losses.size(); ++i)
    {
        if (i != 0) json << ", ";
        json << losses[i];
    }
    json << "],\n"
        << "    \"batch\": " << options.batch << ",\n"
        << "    \"reuse\": " << options.reuse << ",\n"
        << "    \"iterations\": " << options.iterations << ",\n"
        << "    \"warmup\": " << options.warmup << ",\n"
        << "    \"thread_count\": " << options.threads << ",\n"
        << "    \"seed\": " << options.seed << "\n"
        << "  },\n"
        << "  \"resolved\": {\n"
        << "    \"profile\": \"legacy_high_v1\",\n"
        << "    \"field\": \"" << (geometry.parent <= 256 ? "gf8" : "gf16") << "\",\n"
        << "    \"thread_count\": " << resolved_threads << ",\n"
        << "    \"parent_count\": " << geometry.parent << ",\n"
        << "    \"padded_side\": " << geometry.padded_r << "\n"
        << "  },\n"
        << "  \"correctness\": {\"round_trip\": true},\n"
        << "  \"workload_digests\": {\n"
        << "    \"algorithm\": \"fnv1a64\",\n"
        << "    \"original_data\": \"" << original_digest.Hex() << "\",\n"
        << "    \"transmitted_parity\": \"" << parity_digest.Hex() << "\",\n"
        << "    \"recovered_originals\": \"" << recovered_digest.Hex() << "\"\n"
        << "  },\n"
        << "  \"memory\": {\n"
        << "    \"alignment\": 64,\n"
        << "    \"encode_work_shards_per_stripe\": " << encode_count << ",\n"
        << "    \"decode_work_shards_per_stripe\": " << decode_count << ",\n"
        << "    \"encode_work_bytes_batch\": "
        << CheckedProduct(CheckedProduct(encode_count, options.bytes, "encode work"), batch, "encode work") << ",\n"
        << "    \"decode_work_bytes_batch\": "
        << CheckedProduct(CheckedProduct(decode_count, options.bytes, "decode work"), batch, "decode work") << "\n"
        << "  },\n"
        << "  \"metrics\": {\n"
        << "    \"codec_setup\": null,\n"
        << "    \"decode_timing_includes_setup\": true,\n"
        << "    \"encode_execution\": ";
    WriteSummary(json, encode, encode_input, "input_GB_per_s",
        encode_output, "parity_output_GB_per_s", 4);
    json << ",\n    \"decode_including_setup\": ";
    WriteSummary(json, decode, decode_input, "offered_received_GB_per_s",
        decode_output, "repaired_output_GB_per_s", 4);
    json << ",\n    \"rate_semantics\": "
        << "\"offered_received counts every non-null original and all R supplied parity shards\"\n"
        << "  }\n"
        << "}\n";
    return json.str();
}

} // namespace

int main(int argc, char** argv)
{
    try
    {
        const Options options = ParseOptions(argc, argv);
        const std::string json = Run(options);
        if (options.output == "-")
            std::cout << json;
        else
        {
            std::ofstream output(options.output.c_str(),
                std::ios::out | std::ios::trunc);
            if (!output)
                Fail("cannot open JSON output: " + options.output);
            output << json;
            if (!output)
                Fail("failed writing JSON output: " + options.output);
        }
        return 0;
    }
    catch (const std::exception& error)
    {
        std::cerr << "Leopard main benchmark: " << error.what() << std::endl;
        return 1;
    }
}
