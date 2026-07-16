/*
    Copyright (c) 2017 Christopher A. Taylor.  All rights reserved.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the conditions in the top-level
    Leopard LICENSE are met.

    This default-off experiment calls the public Intel ISA-L erasure-code API.
    ISA-L's BSD-3-Clause notice is retained in this directory's NOTICE file.
*/

#include <isa-l/erasure_code.h>

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
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#if !defined(LEO2_ISAL_COMMIT) || !defined(LEO2_ISAL_TREE) || \
    !defined(LEO2_ISAL_LIBRARY_SHA256) || \
    !defined(LEO2_ISAL_HEADER_SHA256) || \
    !defined(LEO2_ISAL_LICENSE_SHA256) || \
    !defined(LEO2_NASM_EXECUTABLE_SHA256) || \
    !defined(LEO2_NASM_ARCHIVE_SHA256)
#error "The standalone CMake identity checks must define ISA-L provenance"
#endif

namespace {

struct Options
{
    uint32_t k;
    uint32_t r;
    uint32_t losses;
    uint64_t bytes;
    size_t batch;
    size_t reuse;
    size_t iterations;
    size_t warmup;
    uint32_t threads;
    uint64_t seed;
    std::string profile;
    std::string oracle;
    std::string output;

    Options()
        : k(240), r(16), losses(1), bytes(65536), batch(1), reuse(8),
          iterations(9), warmup(2), threads(1), seed(1), profile("high"),
          oracle("full"), output("-")
    {}
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

class AlignedBuffer
{
public:
    AlignedBuffer() : pointer_(NULL), bytes_(0) {}
    explicit AlignedBuffer(size_t bytes) : pointer_(NULL), bytes_(0)
    {
        Reset(bytes);
    }
    ~AlignedBuffer() { free(pointer_); }

    void Reset(size_t bytes)
    {
        free(pointer_);
        pointer_ = NULL;
        bytes_ = 0;
        if (bytes == 0)
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

struct Summary
{
    std::vector<double> samples_us;
    double median_us;
    double mad_us;
    double minimum_us;
    double maximum_us;
};

struct EncodePlan
{
    std::vector<uint8_t> matrix;
    std::vector<uint8_t> tables;
};

struct DecodePlan
{
    std::vector<uint8_t> selected_rows;
    std::vector<uint8_t> coefficients;
    std::vector<uint8_t> tables;
};

struct Stripe
{
    AlignedBuffer pristine;
    AlignedBuffer fragments;
    AlignedBuffer restored;
    std::vector<uint8_t*> encode_sources;
    std::vector<uint8_t*> encode_outputs;
    std::vector<uint8_t*> decode_sources;
    std::vector<uint8_t*> decode_outputs;
};

static void Fail(const std::string& message)
{
    throw std::runtime_error(message);
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

static uint32_t ParseU32(const std::string& text, const char* option)
{
    const uint64_t value = ParseUnsigned(text, option);
    if (value > std::numeric_limits<uint32_t>::max())
        Fail(std::string(option) + " exceeds uint32");
    return static_cast<uint32_t>(value);
}

static size_t ParseSize(const std::string& text, const char* option)
{
    const uint64_t value = ParseUnsigned(text, option);
    if (value > std::numeric_limits<size_t>::max())
        Fail(std::string(option) + " exceeds size_t");
    return static_cast<size_t>(value);
}

static const char* NeedValue(int argc, char** argv, int& index)
{
    if (++index >= argc)
        Fail(std::string("missing value for ") + argv[index - 1]);
    return argv[index];
}

static Options ParseOptions(int argc, char** argv)
{
    Options options;
    for (int i = 1; i < argc; ++i)
    {
        const std::string argument(argv[i]);
        if (argument == "--k") options.k = ParseU32(NeedValue(argc, argv, i), "--k");
        else if (argument == "--r") options.r = ParseU32(NeedValue(argc, argv, i), "--r");
        else if (argument == "--loss" || argument == "--losses")
            options.losses = ParseU32(NeedValue(argc, argv, i), "--loss");
        else if (argument == "--bytes")
            options.bytes = ParseUnsigned(NeedValue(argc, argv, i), "--bytes");
        else if (argument == "--batch")
            options.batch = ParseSize(NeedValue(argc, argv, i), "--batch");
        else if (argument == "--reuse")
            options.reuse = ParseSize(NeedValue(argc, argv, i), "--reuse");
        else if (argument == "--iterations")
            options.iterations = ParseSize(NeedValue(argc, argv, i), "--iterations");
        else if (argument == "--warmup")
            options.warmup = ParseSize(NeedValue(argc, argv, i), "--warmup");
        else if (argument == "--threads" || argument == "--thread-count")
            options.threads = ParseU32(NeedValue(argc, argv, i), "--threads");
        else if (argument == "--seed")
            options.seed = ParseUnsigned(NeedValue(argc, argv, i), "--seed");
        else if (argument == "--profile")
            options.profile = NeedValue(argc, argv, i);
        else if (argument == "--oracle")
            options.oracle = NeedValue(argc, argv, i);
        else if (argument == "--json" || argument == "--output")
            options.output = NeedValue(argc, argv, i);
        else if (argument == "--help" || argument == "-h")
        {
            std::cout
                << "Usage: " << argv[0]
                << " --k K --r R --bytes B --loss L --batch N --reuse N"
                << " --iterations N --warmup N --threads 1 --seed N"
                << " --profile auto|high|low --oracle full|projection --json PATH\n";
            exit(0);
        }
        else
            Fail("unknown argument: " + argument);
    }
    if (options.k == 0 || options.r == 0 ||
        static_cast<uint64_t>(options.k) + options.r > 256)
        Fail("ISA-L comparison requires positive K,R and K+R <= 256");
    if (options.bytes == 0 || options.bytes > INT32_MAX ||
        options.bytes > std::numeric_limits<size_t>::max())
        Fail("--bytes must be in 1..INT32_MAX and fit size_t");
    if (options.losses > options.k || options.losses > options.r)
        Fail("--loss must be in 0..min(K,R)");
    if (options.batch == 0 || options.reuse == 0 || options.iterations < 3)
        Fail("--batch and --reuse must be positive; --iterations must be at least 3");
    if (options.threads != 1)
        Fail("this bounded adapter supports exactly one execution thread");
    if (options.profile != "auto" && options.profile != "high" &&
        options.profile != "legacy_high_v1" && options.profile != "low" &&
        options.profile != "low_v1")
        Fail("--profile must be auto, high, legacy_high_v1, low, or low_v1");
    if (options.oracle != "full" && options.oracle != "projection")
        Fail("--oracle must be full or projection");
    return options;
}

static size_t CheckedBytes(uint64_t count, uint64_t bytes, const char* what)
{
    if (count != 0 && bytes > std::numeric_limits<size_t>::max() / count)
        Fail(std::string(what) + " overflows size_t");
    return static_cast<size_t>(count * bytes);
}

static uint64_t CheckedProduct(uint64_t left, uint64_t right, const char* what)
{
    if (left != 0 && right > std::numeric_limits<uint64_t>::max() / left)
        Fail(std::string(what) + " overflows uint64");
    return left * right;
}

static std::vector<uint32_t> SelectLosses(const Options& options)
{
    std::vector<uint32_t> order(options.k);
    for (uint32_t i = 0; i < options.k; ++i)
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

static bool IsLost(const std::vector<uint32_t>& losses, uint32_t index)
{
    return std::binary_search(losses.begin(), losses.end(), index);
}

static EncodePlan BuildEncodePlan(uint32_t k, uint32_t r)
{
    EncodePlan plan;
    plan.matrix.resize(static_cast<size_t>(k + r) * k);
    plan.tables.resize(static_cast<size_t>(32) * k * r);
    gf_gen_cauchy1_matrix(&plan.matrix[0], static_cast<int>(k + r),
        static_cast<int>(k));
    for (uint32_t row = 0; row < k; ++row)
        for (uint32_t column = 0; column < k; ++column)
            if (plan.matrix[static_cast<size_t>(row) * k + column] !=
                static_cast<uint8_t>(row == column))
                Fail("gf_gen_cauchy1_matrix did not produce a systematic prefix");
    ec_init_tables(static_cast<int>(k), static_cast<int>(r),
        &plan.matrix[static_cast<size_t>(k) * k], &plan.tables[0]);
    return plan;
}

static uint8_t ScalarGf8Multiply(uint8_t left, uint8_t right)
{
    uint8_t result = 0;
    for (unsigned bit = 0; bit < 8; ++bit)
    {
        if ((right & 1u) != 0)
            result ^= left;
        const bool high = (left & 0x80u) != 0;
        left = static_cast<uint8_t>(left << 1);
        if (high)
            left ^= 0x1d;
        right = static_cast<uint8_t>(right >> 1);
    }
    return result;
}

static uint8_t ScalarGf8Inverse(uint8_t value)
{
    if (value == 0)
        Fail("cannot invert zero in the independent GF(256) oracle");

    // Fermat inversion, value^(256-2), over x^8+x^4+x^3+x^2+1.  This is
    // deliberately independent of ISA-L's gf_inv tables and generator matrix.
    uint8_t result = 1;
    uint8_t power = value;
    unsigned exponent = 254;
    while (exponent != 0)
    {
        if ((exponent & 1u) != 0)
            result = ScalarGf8Multiply(result, power);
        power = ScalarGf8Multiply(power, power);
        exponent >>= 1;
    }
    return result;
}

static uint8_t IndependentCauchyCoefficient(
    uint32_t k,
    uint32_t parity,
    uint32_t source)
{
    const uint8_t row = static_cast<uint8_t>(k + parity);
    const uint8_t column = static_cast<uint8_t>(source);
    return ScalarGf8Inverse(static_cast<uint8_t>(row ^ column));
}

static std::vector<uint8_t> BuildIndependentCauchyCoefficients(
    const Options& options)
{
    std::vector<uint8_t> coefficients(
        static_cast<size_t>(options.k) * options.r);
    for (uint32_t parity = 0; parity < options.r; ++parity)
        for (uint32_t source = 0; source < options.k; ++source)
            coefficients[static_cast<size_t>(parity) * options.k + source] =
                IndependentCauchyCoefficient(options.k, parity, source);
    return coefficients;
}

static void VerifyIndependentCauchyCoefficients(
    const EncodePlan& encode,
    const Options& options,
    const std::vector<uint8_t>& coefficients)
{
    for (uint32_t parity = 0; parity < options.r; ++parity)
        for (uint32_t source = 0; source < options.k; ++source)
        {
            const uint8_t actual = encode.matrix[
                static_cast<size_t>(options.k + parity) * options.k + source];
            const uint8_t expected = coefficients[
                static_cast<size_t>(parity) * options.k + source];
            if (actual != expected)
                Fail("ISA-L Cauchy coefficient differs from the independent formula");
        }
}

static std::vector<size_t> IndependentParityPositions(const Options& options)
{
    const size_t bytes = static_cast<size_t>(options.bytes);
    std::vector<size_t> positions;
    if (options.oracle == "full")
    {
        positions.resize(bytes);
        for (size_t i = 0; i < bytes; ++i)
            positions[i] = i;
        return positions;
    }

    static const size_t boundaries[] = {
        0, 1, 2, 3, 7, 8, 15, 16, 17, 31, 32, 33, 63, 64, 65,
        127, 128, 129, 255, 256, 257, 1023, 1024, 1025,
    };
    for (size_t i = 0; i < sizeof(boundaries) / sizeof(boundaries[0]); ++i)
        if (boundaries[i] < bytes)
            positions.push_back(boundaries[i]);
    if (bytes > 0)
        positions.push_back(bytes - 1);
    std::sort(positions.begin(), positions.end());
    positions.erase(std::unique(positions.begin(), positions.end()), positions.end());

    XorShift64 random(options.seed ^ UINT64_C(0x4341554348594f52));
    const size_t target = std::min<size_t>(bytes, 64);
    while (positions.size() < target)
    {
        const size_t candidate = static_cast<size_t>(random.Next() % bytes);
        if (std::find(positions.begin(), positions.end(), candidate) == positions.end())
            positions.push_back(candidate);
    }
    std::sort(positions.begin(), positions.end());
    positions.erase(std::unique(positions.begin(), positions.end()), positions.end());
    return positions;
}

static void VerifyIndependentEncoding(
    const Stripe& stripe,
    const Options& options,
    const std::vector<size_t>& positions,
    const std::vector<uint8_t>& coefficients)
{
    const size_t bytes = static_cast<size_t>(options.bytes);
    const size_t original_bytes = CheckedBytes(options.k, options.bytes,
        "pristine original data");
    if (memcmp(stripe.pristine.data(), stripe.fragments.data(), original_bytes) != 0)
        Fail("ISA-L encoding modified a systematic source shard");
    for (uint32_t parity = 0; parity < options.r; ++parity)
    {
        const uint8_t* actual = stripe.fragments.data() +
            static_cast<size_t>(options.k + parity) * bytes;
        for (size_t position_i = 0; position_i < positions.size(); ++position_i)
        {
            const size_t byte = positions[position_i];
            uint8_t expected = 0;
            for (uint32_t source = 0; source < options.k; ++source)
                expected ^= ScalarGf8Multiply(
                    coefficients[static_cast<size_t>(parity) * options.k + source],
                    stripe.pristine.data()[static_cast<size_t>(source) * bytes + byte]);
            if (actual[byte] != expected)
                Fail("ISA-L parity differs from the independent scalar Cauchy oracle");
        }
    }
}

static void VerifyIndependentRecovery(
    const Stripe& stripe,
    const Options& options,
    const std::vector<uint32_t>& losses)
{
    const size_t bytes = static_cast<size_t>(options.bytes);
    const size_t original_bytes = CheckedBytes(
        options.k, options.bytes, "post-decode original data");
    if (memcmp(stripe.pristine.data(), stripe.fragments.data(), original_bytes) != 0)
        Fail("ISA-L execution modified a systematic source shard");
    for (size_t loss_i = 0; loss_i < losses.size(); ++loss_i)
    {
        const uint8_t* expected = stripe.pristine.data() +
            static_cast<size_t>(losses[loss_i]) * bytes;
        const uint8_t* actual = stripe.restored.data() + loss_i * bytes;
        if (memcmp(expected, actual, bytes) != 0)
            Fail("ISA-L recovery differs from immutable deterministic source data");
    }
}

static DecodePlan BuildDecodePlan(
    const EncodePlan& encode,
    uint32_t k,
    uint32_t r,
    const std::vector<uint32_t>& losses)
{
    DecodePlan plan;
    plan.selected_rows.reserve(k);
    for (uint32_t row = 0; row < k && plan.selected_rows.size() < k; ++row)
        if (!IsLost(losses, row))
            plan.selected_rows.push_back(static_cast<uint8_t>(row));
    for (uint32_t parity = 0; parity < r && plan.selected_rows.size() < k; ++parity)
        plan.selected_rows.push_back(static_cast<uint8_t>(k + parity));
    if (plan.selected_rows.size() != k)
        Fail("not enough deterministic received rows for decode");

    std::vector<uint8_t> selected(static_cast<size_t>(k) * k);
    std::vector<uint8_t> inverse(static_cast<size_t>(k) * k);
    for (uint32_t row = 0; row < k; ++row)
        memcpy(&selected[static_cast<size_t>(row) * k],
            &encode.matrix[static_cast<size_t>(plan.selected_rows[row]) * k], k);
    if (gf_invert_matrix(&selected[0], &inverse[0], static_cast<int>(k)) < 0)
        Fail("ISA-L selected received matrix is singular");

    plan.coefficients.resize(static_cast<size_t>(losses.size()) * k);
    for (size_t i = 0; i < losses.size(); ++i)
        memcpy(&plan.coefficients[i * k],
            &inverse[static_cast<size_t>(losses[i]) * k], k);
    plan.tables.resize(static_cast<size_t>(32) * k * losses.size());
    if (!losses.empty())
        ec_init_tables(static_cast<int>(k), static_cast<int>(losses.size()),
            &plan.coefficients[0], &plan.tables[0]);
    return plan;
}

static void InitializeStripe(
    Stripe& stripe,
    const Options& options,
    const DecodePlan& decode,
    size_t stripe_index)
{
    const size_t bytes = static_cast<size_t>(options.bytes);
    stripe.pristine.Reset(CheckedBytes(options.k, options.bytes,
        "pristine original storage"));
    stripe.fragments.Reset(CheckedBytes(options.k + options.r, options.bytes,
        "fragment storage"));
    stripe.restored.Reset(CheckedBytes(options.losses, options.bytes,
        "restored storage"));
    stripe.encode_sources.resize(options.k);
    stripe.encode_outputs.resize(options.r);
    stripe.decode_sources.resize(options.k);
    stripe.decode_outputs.resize(options.losses);

    XorShift64 random(options.seed ^
        (UINT64_C(0x9e3779b97f4a7c15) * (stripe_index + 1)));
    const size_t original_bytes = CheckedBytes(options.k, options.bytes,
        "original data");
    for (size_t i = 0; i < original_bytes; ++i)
        stripe.pristine.data()[i] = static_cast<uint8_t>(random.Next() >> 56);
    memcpy(stripe.fragments.data(), stripe.pristine.data(), original_bytes);
    for (uint32_t i = 0; i < options.k; ++i)
        stripe.encode_sources[i] = stripe.fragments.data() + static_cast<size_t>(i) * bytes;
    for (uint32_t i = 0; i < options.r; ++i)
        stripe.encode_outputs[i] = stripe.fragments.data() +
            static_cast<size_t>(options.k + i) * bytes;
    for (uint32_t i = 0; i < options.k; ++i)
        stripe.decode_sources[i] = stripe.fragments.data() +
            static_cast<size_t>(decode.selected_rows[i]) * bytes;
    for (uint32_t i = 0; i < options.losses; ++i)
        stripe.decode_outputs[i] = stripe.restored.data() + static_cast<size_t>(i) * bytes;
}

static double Median(std::vector<double> values)
{
    const size_t middle = values.size() / 2;
    std::nth_element(values.begin(), values.begin() + middle, values.end());
    const double upper = values[middle];
    if (values.size() & 1u)
        return upper;
    std::nth_element(values.begin(), values.begin() + middle - 1,
        values.begin() + middle);
    return (values[middle - 1] + upper) * 0.5;
}

static Summary Summarize(const std::vector<double>& samples)
{
    Summary result;
    result.samples_us = samples;
    result.median_us = Median(samples);
    std::vector<double> deviations(samples.size());
    for (size_t i = 0; i < samples.size(); ++i)
        deviations[i] = std::fabs(samples[i] - result.median_us);
    result.mad_us = Median(deviations);
    result.minimum_us = *std::min_element(samples.begin(), samples.end());
    result.maximum_us = *std::max_element(samples.begin(), samples.end());
    return result;
}

template<class Callable>
static Summary Measure(size_t iterations, size_t inner_calls, const Callable& callable)
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
        const double elapsed = std::chrono::duration_cast<
            std::chrono::duration<double, std::micro> >(end - begin).count();
        samples.push_back(elapsed / static_cast<double>(inner_calls));
    }
    return Summarize(samples);
}

static uint32_t CeilPow2(uint32_t value)
{
    uint32_t result = 1;
    while (result < value)
        result <<= 1;
    return result;
}

static std::string ResolvedProfile(const Options& options)
{
    if (options.profile == "low" || options.profile == "low_v1")
        return "low_v1";
    if (options.profile == "high" || options.profile == "legacy_high_v1")
        return "legacy_high_v1";
    return options.r > options.k ? "low_v1" : "legacy_high_v1";
}

static uint32_t LeopardParent(const Options& options, const std::string& profile)
{
    if (profile == "low_v1")
        return CeilPow2(CeilPow2(options.k) + options.r);
    return CeilPow2(options.k + CeilPow2(options.r));
}

static double Rate(uint64_t bytes, double microseconds)
{
    return bytes == 0 || microseconds <= 0 ? 0.0 :
        static_cast<double>(bytes) / (microseconds * 1000.0);
}

static void WriteRate(std::ostream& output, uint64_t bytes, double microseconds)
{
    if (bytes == 0)
        output << "null";
    else
        output << Rate(bytes, microseconds);
}

static void WriteArray(std::ostream& output, const std::vector<double>& values)
{
    output << '[';
    for (size_t i = 0; i < values.size(); ++i)
    {
        if (i) output << ',';
        output << values[i];
    }
    output << ']';
}

static void WriteSummary(
    std::ostream& output,
    const Summary& summary,
    uint64_t input_bytes,
    const char* input_rate,
    uint64_t output_bytes,
    const char* output_rate,
    bool execution)
{
    output << "{\"median_us" << (execution ? "_per_batch_call" : "") << "\":"
           << summary.median_us
           << ",\"mad_us" << (execution ? "_per_batch_call" : "") << "\":"
           << summary.mad_us
           << ",\"minimum_us" << (execution ? "_per_batch_call" : "") << "\":"
           << summary.minimum_us
           << ",\"maximum_us" << (execution ? "_per_batch_call" : "") << "\":"
           << summary.maximum_us << ",\"samples_us";
    if (execution) output << "_per_batch_call";
    output << "\":";
    WriteArray(output, summary.samples_us);
    if (execution)
    {
        output << ",\"" << input_rate << "\":";
        WriteRate(output, input_bytes, summary.median_us);
        output << ",\"" << output_rate << "\":";
        WriteRate(output, output_bytes, summary.median_us);
    }
    output << '}';
}

static int Run(const Options& options)
{
    const std::vector<uint32_t> losses = SelectLosses(options);
    const EncodePlan encode = BuildEncodePlan(options.k, options.r);
    const std::vector<uint8_t> independent_coefficients =
        BuildIndependentCauchyCoefficients(options);
    VerifyIndependentCauchyCoefficients(encode, options, independent_coefficients);
    const DecodePlan decode = BuildDecodePlan(encode, options.k, options.r, losses);
    const std::vector<size_t> parity_positions = IndependentParityPositions(options);

    std::vector<Stripe*> stripes;
    stripes.reserve(options.batch);
    try
    {
        for (size_t i = 0; i < options.batch; ++i)
        {
            Stripe* stripe = new Stripe;
            stripes.push_back(stripe);
            InitializeStripe(*stripe, options, decode, i);
        }

        const int length = static_cast<int>(options.bytes);
        const int k = static_cast<int>(options.k);
        const int r = static_cast<int>(options.r);
        const int lost = static_cast<int>(options.losses);
        const auto encode_call = [&]() {
            for (size_t i = 0; i < stripes.size(); ++i)
                ec_encode_data(length, k, r, const_cast<uint8_t*>(&encode.tables[0]),
                    &stripes[i]->encode_sources[0], &stripes[i]->encode_outputs[0]);
        };
        const auto decode_call = [&]() {
            if (lost > 0)
                for (size_t i = 0; i < stripes.size(); ++i)
                    ec_encode_data(length, k, lost,
                        const_cast<uint8_t*>(&decode.tables[0]),
                        &stripes[i]->decode_sources[0],
                        &stripes[i]->decode_outputs[0]);
        };

        encode_call();
        for (size_t stripe_i = 0; stripe_i < stripes.size(); ++stripe_i)
            VerifyIndependentEncoding(
                *stripes[stripe_i], options, parity_positions,
                independent_coefficients);
        decode_call();
        for (size_t stripe_i = 0; stripe_i < stripes.size(); ++stripe_i)
            VerifyIndependentRecovery(*stripes[stripe_i], options, losses);

        const Summary codec_setup = Measure(options.iterations, 1, [&]() {
            const EncodePlan temporary = BuildEncodePlan(options.k, options.r);
            if (temporary.tables.empty()) Fail("empty encode setup");
        });
        const Summary plan_setup = Measure(options.iterations, 1, [&]() {
            const DecodePlan temporary = BuildDecodePlan(
                encode, options.k, options.r, losses);
            if (options.losses > 0 && temporary.tables.empty())
                Fail("empty decode setup");
        });
        for (size_t i = 0; i < options.warmup; ++i)
        {
            encode_call();
            decode_call();
        }
        const Summary encode_execution = Measure(
            options.iterations, options.reuse, encode_call);
        const Summary decode_execution = Measure(
            options.iterations, options.reuse, decode_call);
        for (size_t stripe_i = 0; stripe_i < stripes.size(); ++stripe_i)
        {
            VerifyIndependentEncoding(
                *stripes[stripe_i], options, parity_positions,
                independent_coefficients);
            VerifyIndependentRecovery(*stripes[stripe_i], options, losses);
        }

        const uint64_t batch = static_cast<uint64_t>(options.batch);
        const uint64_t encode_input = CheckedProduct(
            CheckedProduct(options.k, options.bytes, "encode input"), batch,
            "encode input");
        const uint64_t encode_output = CheckedProduct(
            CheckedProduct(options.r, options.bytes, "encode output"), batch,
            "encode output");
        const uint64_t offered_shards = options.k - options.losses + options.r;
        const uint64_t decode_offered = CheckedProduct(
            CheckedProduct(offered_shards, options.bytes, "decode offered"), batch,
            "decode offered");
        const uint64_t decode_selected = CheckedProduct(
            CheckedProduct(options.k, options.bytes, "decode selected"), batch,
            "decode selected");
        const uint64_t decode_output = CheckedProduct(
            CheckedProduct(options.losses, options.bytes, "decode output"), batch,
            "decode output");
        const uint64_t coefficient_checks = CheckedProduct(
            options.k, options.r, "independent coefficient checks");
        const uint64_t parity_checked = CheckedProduct(
            CheckedProduct(static_cast<uint64_t>(parity_positions.size()), options.r,
                "independent parity checked bytes"), batch,
            "independent parity checked bytes");
        const uint64_t parity_total = CheckedProduct(
            CheckedProduct(options.bytes, options.r,
                "independent parity total bytes"), batch,
            "independent parity total bytes");
        const std::string profile = ResolvedProfile(options);
        const uint32_t parent = LeopardParent(options, profile);
        const char* leopard_field = parent <= 256 ? "gf8" : "gf16";

        std::ostringstream json;
        json << std::fixed << std::setprecision(9);
        json << "{\n"
             << "  \"schema\":\"leopard2-isal-benchmark/v2\",\n"
             << "  \"provider\":{\"name\":\"Intel ISA-L\",\"source_commit\":\""
             << LEO2_ISAL_COMMIT << "\",\"source_tree\":\""
             << LEO2_ISAL_TREE << "\",\"library_sha256\":\""
             << LEO2_ISAL_LIBRARY_SHA256
             << "\",\"header_sha256\":\""
             << LEO2_ISAL_HEADER_SHA256
             << "\",\"license\":\"BSD-3-Clause\",\"license_sha256\":\""
             << LEO2_ISAL_LICENSE_SHA256
             << "\",\"nasm_executable_sha256\":\""
             << LEO2_NASM_EXECUTABLE_SHA256
             << "\",\"nasm_archive_sha256\":\""
             << LEO2_NASM_ARCHIVE_SHA256 << "\",\"field\":\"gf8\","
             << "\"generator\":\"gf_gen_cauchy1_matrix\",\"wire_compatible\":false},\n"
             << "  \"parameters\":{\"K\":" << options.k << ",\"R\":" << options.r
             << ",\"requested_profile\":\"" << profile
             << "\",\"requested_field\":\"auto\",\"requested_backend\":\"auto\""
             << ",\"shard_bytes\":" << options.bytes
             << ",\"loss_count\":" << options.losses
             << ",\"missing_original_indices\":[";
        for (size_t i = 0; i < losses.size(); ++i)
        {
            if (i) json << ',';
            json << losses[i];
        }
        json << "],\"batch\":" << options.batch << ",\"reuse\":" << options.reuse
             << ",\"iterations\":" << options.iterations << ",\"warmup\":"
             << options.warmup << ",\"thread_count\":1,\"seed\":" << options.seed
             << "},\n"
             << "  \"comparison_identity\":{\"leopard2_profile\":\"" << profile
             << "\",\"leopard2_parent\":" << parent
             << ",\"leopard2_field\":\"" << leopard_field
             << "\",\"isa_l_field_advantage_from_padding\":"
             << ((parent > 256 &&
                  static_cast<uint64_t>(options.k) + options.r <= 256)
                    ? "true" : "false")
             << ",\"scope\":\"public payload and repaired-output throughput only; parity bytes and generator matrices differ\"},\n"
             << "  \"correctness\":{\"direct_source_round_trip\":true,"
             << "\"systematic_generator_prefix\":true,"
             << "\"systematic_sources_immutable\":true,"
             << "\"independent_scalar_cauchy_coefficients_checked\":"
             << coefficient_checks << ','
             << "\"independent_scalar_parity_mode\":\"" << options.oracle << "\","
             << "\"independent_scalar_parity_checked_bytes_per_validation\":"
             << parity_checked
             << ",\"independent_scalar_parity_total_bytes_per_validation\":"
             << parity_total
             << ",\"independent_scalar_parity_validation_passes\":2},\n"
             << "  \"memory\":{\"alignment_bytes\":64,\"direct_application_buffers\":true,"
             << "\"staging_copy_bytes_per_execution\":0,\"encode_input_bytes_per_batch_call\":"
             << encode_input << ",\"encode_output_bytes_per_batch_call\":" << encode_output
             << ",\"decode_offered_bytes_per_batch_call\":" << decode_offered
             << ",\"decode_selected_bytes_per_batch_call\":" << decode_selected
             << ",\"decode_output_bytes_per_batch_call\":" << decode_output
             << ",\"matrix_bytes\":" << encode.matrix.size()
             << ",\"encode_table_bytes\":" << encode.tables.size()
             << ",\"decode_table_bytes\":" << decode.tables.size() << "},\n"
             << "  \"metrics\":{\"codec_setup\":";
        WriteSummary(json, codec_setup, 0, "", 0, "", false);
        json << ",\"encode_execution\":";
        WriteSummary(json, encode_execution, encode_input, "input_GB_per_s",
            encode_output, "parity_output_GB_per_s", true);
        json << ",\"decode_plan_setup\":";
        WriteSummary(json, plan_setup, 0, "", 0, "", false);
        json << ",\"decode_execution\":";
        WriteSummary(json, decode_execution, decode_offered,
            "offered_received_GB_per_s", decode_output,
            "repaired_output_GB_per_s", true);
        const double amortized = decode_execution.median_us +
            plan_setup.median_us / static_cast<double>(options.reuse);
        json << ",\"decode_amortized_at_reuse\":{\"reuse_count\":" << options.reuse
             << ",\"derived_median_us_per_batch_call\":" << amortized
             << ",\"offered_received_GB_per_s\":" << Rate(decode_offered, amortized)
             << ",\"repaired_output_GB_per_s\":";
        WriteRate(json, decode_output, amortized);
        json
             << "},\"rate_semantics\":\"offered_received counts all non-erased public shards; ISA-L reads the recorded deterministic K-row subset\"}\n"
             << "}\n";

        if (options.output == "-")
            std::cout << json.str();
        else
        {
            std::ofstream output(options.output.c_str(), std::ios::binary);
            if (!output) Fail("cannot open output: " + options.output);
            output << json.str();
            if (!output) Fail("cannot write output: " + options.output);
        }
    }
    catch (...)
    {
        for (size_t i = 0; i < stripes.size(); ++i) delete stripes[i];
        throw;
    }
    for (size_t i = 0; i < stripes.size(); ++i) delete stripes[i];
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
        std::cerr << "leopard2_isal_benchmark: " << error.what() << '\n';
        return 1;
    }
}
