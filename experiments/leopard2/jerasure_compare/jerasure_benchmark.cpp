/*
    Copyright (c) 2017 Christopher A. Taylor.  All rights reserved.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the conditions in the top-level
    Leopard LICENSE are met.

    This default-off executable calls the public Jerasure 2.0 API.  The exact
    Jerasure and GF-Complete BSD notices are retained in this directory.
*/

#include <jerasure.h>
#include <reed_sol.h>

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

#if !defined(LEO2_JERASURE_COMMIT) || !defined(LEO2_JERASURE_TREE) || \
    !defined(LEO2_JERASURE_HEADER_SHA256) || \
    !defined(LEO2_REED_SOL_HEADER_SHA256) || \
    !defined(LEO2_JERASURE_LICENSE_SHA256) || \
    !defined(LEO2_GF_COMPLETE_COMMIT) || !defined(LEO2_GF_COMPLETE_TREE) || \
    !defined(LEO2_GF_COMPLETE_HEADER_SHA256) || \
    !defined(LEO2_GF_COMPLETE_LICENSE_SHA256) || \
    !defined(LEO2_GF_COMPLETE_SIMD_FLAGS)
#error "The standalone CMake identity checks must define provider provenance"
#endif

namespace {

static_assert(sizeof(int) == 4,
    "Jerasure comparison evidence defines matrix storage as 32-bit int");

static const uint32_t kMaximumOriginals = 4096;
static const uint32_t kMaximumRecovery = 4096;
static const uint64_t kMaximumMatrixCoefficients = UINT64_C(1) << 24;
static const uint64_t kMaximumApplicationBytes = UINT64_C(8) << 30;
static const uint32_t kMaximumLeopardParent = 65536;

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
    explicit AlignedBuffer(size_t bytes) : pointer_(NULL), bytes_(0) { Reset(bytes); }
    ~AlignedBuffer() { free(pointer_); }

    void Reset(size_t bytes)
    {
        free(pointer_);
        pointer_ = NULL;
        bytes_ = 0;
        if (bytes == 0) return;
        if (posix_memalign(&pointer_, 64, bytes) != 0) throw std::bad_alloc();
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
    std::vector<int> matrix;
};

struct DecodePlan
{
    std::vector<int> decoding_matrix;
    std::vector<int> selected_ids;
    std::vector<int> missing_rows;
};

struct Stripe
{
    AlignedBuffer pristine;
    AlignedBuffer fragments;
    AlignedBuffer restored;
    std::vector<char*> encode_data;
    std::vector<char*> decode_data;
    std::vector<char*> coding;
};

static void Fail(const std::string& message) { throw std::runtime_error(message); }

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
    if (++index >= argc) Fail(std::string("missing value for ") + argv[index - 1]);
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
        else if (argument == "--profile") options.profile = NeedValue(argc, argv, i);
        else if (argument == "--oracle") options.oracle = NeedValue(argc, argv, i);
        else if (argument == "--json" || argument == "--output")
            options.output = NeedValue(argc, argv, i);
        else if (argument == "--help" || argument == "-h")
        {
            std::cout << "Usage: " << argv[0]
                << " --k K --r R --bytes B --loss L --batch N --reuse N"
                << " --iterations N --warmup N --threads 1 --seed N"
                << " --profile auto|high|low --oracle full|projection --json PATH\n";
            exit(0);
        }
        else Fail("unknown argument: " + argument);
    }
    const uint64_t total = static_cast<uint64_t>(options.k) + options.r;
    if (options.k == 0 || options.r == 0 || total > 65536)
        Fail("Jerasure comparison requires positive K,R and K+R <= 65536");
    if (options.k > kMaximumOriginals || options.r > kMaximumRecovery)
        Fail("bounded Jerasure comparison requires K,R <= 4096");
    if (static_cast<uint64_t>(options.k) * options.r >
            kMaximumMatrixCoefficients ||
        (options.losses != 0 && static_cast<uint64_t>(options.k) * options.k >
            kMaximumMatrixCoefficients))
        Fail("bounded Jerasure comparison matrix domain exceeded");
    if (options.bytes == 0 || options.bytes > INT32_MAX ||
        options.bytes > std::numeric_limits<size_t>::max())
        Fail("--bytes must be in 1..INT32_MAX and fit size_t");
    if ((options.bytes & 7u) != 0)
        Fail("--bytes must satisfy Jerasure's deterministic 8-byte region contract");
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
    const uint64_t slots_per_stripe =
        static_cast<uint64_t>(options.k) * 2 + options.r + options.losses;
    if (slots_per_stripe > kMaximumApplicationBytes / options.batch ||
        slots_per_stripe * options.batch > kMaximumApplicationBytes / options.bytes)
        Fail("bounded Jerasure comparison application-buffer domain exceeded");
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
    for (uint32_t i = 0; i < options.k; ++i) order[i] = i;
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

static int FieldWidth(const Options& options)
{
    return static_cast<uint64_t>(options.k) + options.r <= 256 ? 8 : 16;
}

static uint32_t PrimitivePolynomial(int width)
{
    return width == 8 ? UINT32_C(0x11d) : UINT32_C(0x1100b);
}

static uint32_t ScalarMultiply(uint32_t left, uint32_t right, int width)
{
    const uint32_t high = UINT32_C(1) << (width - 1);
    const uint32_t mask = (UINT32_C(1) << width) - 1;
    const uint32_t reduction = PrimitivePolynomial(width) & mask;
    uint32_t result = 0;
    for (int bit = 0; bit < width; ++bit)
    {
        if (right & 1u) result ^= left;
        const bool carry = (left & high) != 0;
        left = (left << 1) & mask;
        if (carry) left ^= reduction;
        right >>= 1;
    }
    return result & mask;
}

static uint32_t ScalarInverse(uint32_t value, int width)
{
    if (value == 0) Fail("independent field oracle cannot invert zero");
    uint32_t result = 1;
    uint32_t power = value;
    uint32_t exponent = (UINT32_C(1) << width) - 2;
    while (exponent != 0)
    {
        if (exponent & 1u) result = ScalarMultiply(result, power, width);
        power = ScalarMultiply(power, power, width);
        exponent >>= 1;
    }
    return result;
}

static std::vector<uint32_t> InvertMatrix(
    const std::vector<uint32_t>& input, uint32_t rows, int width)
{
    std::vector<uint32_t> work(input);
    std::vector<uint32_t> inverse(static_cast<size_t>(rows) * rows, 0);
    for (uint32_t i = 0; i < rows; ++i)
        inverse[static_cast<size_t>(i) * rows + i] = 1;
    for (uint32_t column = 0; column < rows; ++column)
    {
        uint32_t pivot = column;
        while (pivot < rows && work[static_cast<size_t>(pivot) * rows + column] == 0)
            ++pivot;
        if (pivot == rows) Fail("independent Vandermonde prefix is singular");
        if (pivot != column)
            for (uint32_t j = 0; j < rows; ++j)
            {
                std::swap(work[static_cast<size_t>(column) * rows + j],
                          work[static_cast<size_t>(pivot) * rows + j]);
                std::swap(inverse[static_cast<size_t>(column) * rows + j],
                          inverse[static_cast<size_t>(pivot) * rows + j]);
            }
        const uint32_t scale = ScalarInverse(
            work[static_cast<size_t>(column) * rows + column], width);
        for (uint32_t j = 0; j < rows; ++j)
        {
            work[static_cast<size_t>(column) * rows + j] = ScalarMultiply(
                work[static_cast<size_t>(column) * rows + j], scale, width);
            inverse[static_cast<size_t>(column) * rows + j] = ScalarMultiply(
                inverse[static_cast<size_t>(column) * rows + j], scale, width);
        }
        for (uint32_t row = 0; row < rows; ++row)
        {
            if (row == column) continue;
            const uint32_t factor = work[static_cast<size_t>(row) * rows + column];
            if (factor == 0) continue;
            for (uint32_t j = 0; j < rows; ++j)
            {
                work[static_cast<size_t>(row) * rows + j] ^= ScalarMultiply(
                    factor, work[static_cast<size_t>(column) * rows + j], width);
                inverse[static_cast<size_t>(row) * rows + j] ^= ScalarMultiply(
                    factor, inverse[static_cast<size_t>(column) * rows + j], width);
            }
        }
    }
    return inverse;
}

static std::vector<uint32_t> ExtendedVandermondeRow(
    uint32_t row, uint32_t rows, uint32_t columns, int width)
{
    std::vector<uint32_t> result(columns, 0);
    if (row == 0) result[0] = 1;
    else if (row + 1 == rows) result[columns - 1] = 1;
    else
    {
        uint32_t power = 1;
        for (uint32_t column = 0; column < columns; ++column)
        {
            result[column] = power;
            power = ScalarMultiply(power, row, width);
        }
    }
    return result;
}

static std::vector<uint32_t> BuildIndependentCodingMatrix(
    uint32_t k, uint32_t r, int width)
{
    const uint32_t rows = k + r;
    std::vector<uint32_t> prefix(static_cast<size_t>(k) * k);
    for (uint32_t row = 0; row < k; ++row)
    {
        const std::vector<uint32_t> values =
            ExtendedVandermondeRow(row, rows, k, width);
        std::copy(values.begin(), values.end(), prefix.begin() + static_cast<size_t>(row) * k);
    }
    const std::vector<uint32_t> inverse = InvertMatrix(prefix, k, width);
    std::vector<uint32_t> result(static_cast<size_t>(r) * k, 0);
    for (uint32_t parity = 0; parity < r; ++parity)
    {
        const std::vector<uint32_t> values =
            ExtendedVandermondeRow(k + parity, rows, k, width);
        for (uint32_t output = 0; output < k; ++output)
            for (uint32_t term = 0; term < k; ++term)
                result[static_cast<size_t>(parity) * k + output] ^=
                    ScalarMultiply(values[term],
                        inverse[static_cast<size_t>(term) * k + output], width);
    }
    // Jerasure's distribution-matrix routine applies two documented GRS row
    // normalizations after systematizing the Vandermonde prefix: parity row
    // zero becomes all ones, then column zero becomes one in every later
    // parity row.  Reproduce those algebraic scalings independently rather
    // than using Jerasure/GF-Complete arithmetic as the oracle.
    for (uint32_t column = 0; column < k; ++column)
    {
        const uint32_t scale = ScalarInverse(result[column], width);
        for (uint32_t parity = 0; parity < r; ++parity)
            result[static_cast<size_t>(parity) * k + column] = ScalarMultiply(
                result[static_cast<size_t>(parity) * k + column], scale, width);
    }
    for (uint32_t parity = 1; parity < r; ++parity)
    {
        const size_t row = static_cast<size_t>(parity) * k;
        const uint32_t scale = ScalarInverse(result[row], width);
        for (uint32_t column = 0; column < k; ++column)
            result[row + column] = ScalarMultiply(result[row + column], scale, width);
    }
    return result;
}

static EncodePlan BuildEncodePlan(uint32_t k, uint32_t r, int width)
{
    int* raw = reed_sol_vandermonde_coding_matrix(
        static_cast<int>(k), static_cast<int>(r), width);
    if (raw == NULL) Fail("reed_sol_vandermonde_coding_matrix failed");
    EncodePlan plan;
    plan.matrix.assign(raw, raw + static_cast<size_t>(k) * r);
    free(raw);
    return plan;
}

static void VerifyIndependentMatrix(
    const EncodePlan& plan, const std::vector<uint32_t>& independent)
{
    if (plan.matrix.size() != independent.size()) Fail("coding-matrix size mismatch");
    for (size_t i = 0; i < independent.size(); ++i)
        if (static_cast<uint32_t>(plan.matrix[i]) != independent[i])
            Fail("Jerasure coding matrix differs from independent systematic Vandermonde algebra");
}

static DecodePlan BuildDecodePlan(
    const EncodePlan& encode, uint32_t k, uint32_t r, int width,
    const std::vector<uint32_t>& losses)
{
    DecodePlan plan;
    if (losses.empty()) return plan;
    plan.decoding_matrix.resize(static_cast<size_t>(k) * k);
    plan.selected_ids.resize(k);
    plan.missing_rows.reserve(losses.size());
    std::vector<int> erased(static_cast<size_t>(k) + r, 0);
    for (size_t i = 0; i < losses.size(); ++i) erased[losses[i]] = 1;
    std::vector<int> mutable_matrix(encode.matrix);
    if (jerasure_make_decoding_matrix(
            static_cast<int>(k), static_cast<int>(r), width,
            &mutable_matrix[0], &erased[0], &plan.decoding_matrix[0],
            &plan.selected_ids[0]) != 0)
        Fail("jerasure_make_decoding_matrix failed");
    for (size_t i = 0; i < losses.size(); ++i)
        plan.missing_rows.push_back(static_cast<int>(losses[i]));
    return plan;
}

static void InitializeStripe(
    Stripe& stripe, const Options& options, const DecodePlan& decode,
    const std::vector<uint32_t>& losses, size_t stripe_index)
{
    const size_t bytes = static_cast<size_t>(options.bytes);
    stripe.pristine.Reset(CheckedBytes(options.k, options.bytes, "pristine originals"));
    stripe.fragments.Reset(CheckedBytes(options.k + options.r, options.bytes, "fragments"));
    stripe.restored.Reset(CheckedBytes(options.losses, options.bytes, "restored originals"));
    stripe.encode_data.resize(options.k);
    stripe.decode_data.resize(options.k);
    stripe.coding.resize(options.r);

    XorShift64 random(options.seed ^
        (UINT64_C(0x9e3779b97f4a7c15) * (stripe_index + 1)));
    for (size_t i = 0; i < stripe.pristine.size(); ++i)
        stripe.pristine.data()[i] = static_cast<uint8_t>(random.Next() >> 56);
    memcpy(stripe.fragments.data(), stripe.pristine.data(), stripe.pristine.size());
    for (uint32_t i = 0; i < options.k; ++i)
    {
        stripe.encode_data[i] = reinterpret_cast<char*>(
            stripe.fragments.data() + static_cast<size_t>(i) * bytes);
        stripe.decode_data[i] = stripe.encode_data[i];
    }
    for (uint32_t i = 0; i < options.r; ++i)
        stripe.coding[i] = reinterpret_cast<char*>(
            stripe.fragments.data() + static_cast<size_t>(options.k + i) * bytes);
    for (size_t i = 0; i < losses.size(); ++i)
        stripe.decode_data[losses[i]] = reinterpret_cast<char*>(
            stripe.restored.data() + i * bytes);
    if (!losses.empty() && decode.selected_ids.empty()) Fail("decode source IDs missing");
}

static uint32_t ReadSymbol(const uint8_t* data, size_t symbol, int width)
{
    if (width == 8) return data[symbol];
    uint16_t value;
    memcpy(&value, data + symbol * 2, sizeof(value));
    return value;
}

static std::vector<size_t> IndependentSymbolPositions(const Options& options, int width)
{
    const size_t symbols = static_cast<size_t>(options.bytes) / (width / 8);
    std::vector<size_t> positions;
    if (options.oracle == "full")
    {
        positions.resize(symbols);
        for (size_t i = 0; i < symbols; ++i) positions[i] = i;
        return positions;
    }
    static const size_t boundaries[] =
        {0, 1, 2, 3, 7, 8, 15, 16, 17, 31, 32, 33, 63, 64, 65,
         127, 128, 129, 255, 256, 257, 511, 512, 513};
    for (size_t i = 0; i < sizeof(boundaries) / sizeof(boundaries[0]); ++i)
        if (boundaries[i] < symbols) positions.push_back(boundaries[i]);
    if (symbols != 0) positions.push_back(symbols - 1);
    std::sort(positions.begin(), positions.end());
    positions.erase(std::unique(positions.begin(), positions.end()), positions.end());
    XorShift64 random(options.seed ^ UINT64_C(0x4a45524153555245));
    const size_t target = std::min<size_t>(symbols, 64 / (width / 8));
    while (positions.size() < target)
    {
        const size_t candidate = static_cast<size_t>(random.Next() % symbols);
        if (std::find(positions.begin(), positions.end(), candidate) == positions.end())
            positions.push_back(candidate);
    }
    std::sort(positions.begin(), positions.end());
    return positions;
}

static void VerifyIndependentEncoding(
    const Stripe& stripe, const Options& options, int width,
    const std::vector<size_t>& positions,
    const std::vector<uint32_t>& coefficients)
{
    const size_t bytes = static_cast<size_t>(options.bytes);
    if (memcmp(stripe.pristine.data(), stripe.fragments.data(), stripe.pristine.size()) != 0)
        Fail("Jerasure encoding modified a systematic source shard");
    for (uint32_t parity = 0; parity < options.r; ++parity)
    {
        const uint8_t* actual = stripe.fragments.data() +
            static_cast<size_t>(options.k + parity) * bytes;
        for (size_t position_i = 0; position_i < positions.size(); ++position_i)
        {
            const size_t symbol = positions[position_i];
            uint32_t expected = 0;
            for (uint32_t source = 0; source < options.k; ++source)
                expected ^= ScalarMultiply(
                    coefficients[static_cast<size_t>(parity) * options.k + source],
                    ReadSymbol(stripe.pristine.data() + static_cast<size_t>(source) * bytes,
                               symbol, width), width);
            if (ReadSymbol(actual, symbol, width) != expected)
                Fail("Jerasure parity differs from independent scalar Vandermonde oracle");
        }
    }
}

static void VerifyRecovery(
    const Stripe& stripe, const Options& options,
    const std::vector<uint32_t>& losses)
{
    const size_t bytes = static_cast<size_t>(options.bytes);
    if (memcmp(stripe.pristine.data(), stripe.fragments.data(), stripe.pristine.size()) != 0)
        Fail("Jerasure execution modified application source storage");
    for (size_t i = 0; i < losses.size(); ++i)
        if (memcmp(stripe.pristine.data() + static_cast<size_t>(losses[i]) * bytes,
                   stripe.restored.data() + i * bytes, bytes) != 0)
            Fail("Jerasure recovery differs from immutable deterministic source data");
}

static double Median(std::vector<double> values)
{
    const size_t middle = values.size() / 2;
    std::nth_element(values.begin(), values.begin() + middle, values.end());
    const double upper = values[middle];
    if (values.size() & 1u) return upper;
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
        for (size_t call = 0; call < inner_calls; ++call) callable();
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
    while (result < value) result <<= 1;
    return result;
}

static std::string ResolvedProfile(const Options& options)
{
    if (options.profile == "low" || options.profile == "low_v1") return "low_v1";
    if (options.profile == "high" || options.profile == "legacy_high_v1")
        return "legacy_high_v1";
    return options.r > options.k ? "low_v1" : "legacy_high_v1";
}

static uint32_t LeopardParent(const Options& options, const std::string& profile)
{
    if (profile == "low_v1") return CeilPow2(CeilPow2(options.k) + options.r);
    return CeilPow2(options.k + CeilPow2(options.r));
}

static double Rate(uint64_t bytes, double microseconds)
{
    return bytes == 0 || microseconds <= 0 ? 0.0 :
        static_cast<double>(bytes) / (microseconds * 1000.0);
}

static void WriteRate(std::ostream& output, uint64_t bytes, double microseconds)
{
    if (bytes == 0) output << "null";
    else output << Rate(bytes, microseconds);
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
    std::ostream& output, const Summary& summary, uint64_t input_bytes,
    const char* input_rate, uint64_t output_bytes, const char* output_rate,
    bool execution)
{
    const char* suffix = execution ? "_per_batch_call" : "";
    output << "{\"median_us" << suffix << "\":" << summary.median_us
           << ",\"mad_us" << suffix << "\":" << summary.mad_us
           << ",\"minimum_us" << suffix << "\":" << summary.minimum_us
           << ",\"maximum_us" << suffix << "\":" << summary.maximum_us
           << ",\"samples_us" << suffix << "\":";
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

static uint64_t Fnv1a64Update(uint64_t digest, const void* data, size_t bytes)
{
    static const uint64_t prime = UINT64_C(1099511628211);
    const uint8_t* input = static_cast<const uint8_t*>(data);
    for (size_t i = 0; i < bytes; ++i)
    {
        digest ^= input[i];
        digest *= prime;
    }
    return digest;
}

static std::string HexU64(uint64_t value)
{
    std::ostringstream output;
    output << std::hex << std::nouppercase << std::setw(16)
           << std::setfill('0') << value;
    return output.str();
}

static int Run(const Options& options)
{
    const int width = FieldWidth(options);
    const std::vector<uint32_t> losses = SelectLosses(options);
    const std::string profile = ResolvedProfile(options);
    const uint32_t parent = LeopardParent(options, profile);
    if (parent > kMaximumLeopardParent)
        Fail("paired Leopard2 parent exceeds 65536 coordinates");
    const EncodePlan encode = BuildEncodePlan(options.k, options.r, width);
    const std::vector<uint32_t> independent =
        BuildIndependentCodingMatrix(options.k, options.r, width);
    VerifyIndependentMatrix(encode, independent);
    const DecodePlan decode = BuildDecodePlan(
        encode, options.k, options.r, width, losses);
    const std::vector<size_t> positions = IndependentSymbolPositions(options, width);

    std::vector<Stripe*> stripes;
    stripes.reserve(options.batch);
    try
    {
        for (size_t i = 0; i < options.batch; ++i)
        {
            Stripe* stripe = new Stripe;
            stripes.push_back(stripe);
            InitializeStripe(*stripe, options, decode, losses, i);
        }
        const int k = static_cast<int>(options.k);
        const int r = static_cast<int>(options.r);
        const int bytes = static_cast<int>(options.bytes);
        const auto encode_call = [&]() {
            for (size_t i = 0; i < stripes.size(); ++i)
                jerasure_matrix_encode(k, r, width,
                    const_cast<int*>(&encode.matrix[0]), &stripes[i]->encode_data[0],
                    &stripes[i]->coding[0], bytes);
        };
        const auto decode_call = [&]() {
            for (size_t stripe_i = 0; stripe_i < stripes.size(); ++stripe_i)
                for (size_t loss_i = 0; loss_i < losses.size(); ++loss_i)
                    jerasure_matrix_dotprod(k, width,
                        const_cast<int*>(&decode.decoding_matrix[
                            static_cast<size_t>(losses[loss_i]) * options.k]),
                        const_cast<int*>(&decode.selected_ids[0]),
                        static_cast<int>(losses[loss_i]),
                        &stripes[stripe_i]->decode_data[0],
                        &stripes[stripe_i]->coding[0], bytes);
        };

        encode_call();
        for (size_t i = 0; i < stripes.size(); ++i)
            VerifyIndependentEncoding(*stripes[i], options, width, positions, independent);
        decode_call();
        for (size_t i = 0; i < stripes.size(); ++i)
            VerifyRecovery(*stripes[i], options, losses);

        static const uint64_t offset = UINT64_C(14695981039346656037);
        uint64_t original_digest = offset;
        uint64_t parity_digest = offset;
        uint64_t recovered_digest = offset;
        const size_t shard_bytes = static_cast<size_t>(options.bytes);
        for (size_t i = 0; i < stripes.size(); ++i)
        {
            original_digest = Fnv1a64Update(original_digest,
                stripes[i]->pristine.data(), stripes[i]->pristine.size());
            parity_digest = Fnv1a64Update(parity_digest,
                stripes[i]->fragments.data() + static_cast<size_t>(options.k) * shard_bytes,
                CheckedBytes(options.r, options.bytes, "parity digest"));
            for (size_t loss_i = 0; loss_i < losses.size(); ++loss_i)
                recovered_digest = Fnv1a64Update(recovered_digest,
                    stripes[i]->restored.data() + loss_i * shard_bytes, shard_bytes);
        }

        const Summary codec_setup = Measure(options.iterations, 1, [&]() {
            const EncodePlan temporary = BuildEncodePlan(options.k, options.r, width);
            if (temporary.matrix.empty()) Fail("empty Jerasure encode setup");
        });
        const Summary plan_setup = Measure(options.iterations, 1, [&]() {
            const DecodePlan temporary = BuildDecodePlan(
                encode, options.k, options.r, width, losses);
            if (!losses.empty() && temporary.selected_ids.empty())
                Fail("empty Jerasure decode setup");
        });
        for (size_t i = 0; i < options.warmup; ++i)
        {
            encode_call();
            decode_call();
        }
        const Summary encode_execution = Measure(options.iterations, options.reuse, encode_call);
        const Summary decode_execution = Measure(options.iterations, options.reuse, decode_call);
        for (size_t i = 0; i < stripes.size(); ++i)
        {
            VerifyIndependentEncoding(*stripes[i], options, width, positions, independent);
            VerifyRecovery(*stripes[i], options, losses);
        }

        const uint64_t batch = options.batch;
        const uint64_t encode_input = CheckedProduct(
            CheckedProduct(options.k, options.bytes, "encode input"), batch, "encode input");
        const uint64_t encode_output = CheckedProduct(
            CheckedProduct(options.r, options.bytes, "encode output"), batch, "encode output");
        const uint64_t offered = losses.empty() ? 0 : CheckedProduct(
            CheckedProduct(options.k - options.losses + options.r, options.bytes,
                           "decode offered"), batch, "decode offered");
        const uint64_t selected = losses.empty() ? 0 : CheckedProduct(
            CheckedProduct(options.k, options.bytes, "decode selected"), batch,
            "decode selected");
        const uint64_t decode_output = CheckedProduct(
            CheckedProduct(options.losses, options.bytes, "decode output"), batch,
            "decode output");
        const uint64_t checked_parity = CheckedProduct(
            CheckedProduct(positions.size() * static_cast<size_t>(width / 8),
                           options.r, "checked parity"), batch, "checked parity");
        const uint64_t total_parity = CheckedProduct(
            CheckedProduct(options.bytes, options.r, "total parity"), batch,
            "total parity");
        const char* leopard_field = parent <= 256 ? "gf8" : "gf16";
        const char* provider_field = width == 8 ? "gf8" : "gf16";

        std::ostringstream json;
        json << std::fixed << std::setprecision(9);
        json << "{\n  \"schema\":\"leopard2-jerasure-benchmark/v1\",\n"
             << "  \"provider\":{\"name\":\"Jerasure 2.0\",\"source_commit\":\""
             << LEO2_JERASURE_COMMIT << "\",\"source_tree\":\""
             << LEO2_JERASURE_TREE << "\",\"header_sha256\":\""
             << LEO2_JERASURE_HEADER_SHA256 << "\",\"reed_sol_header_sha256\":\""
             << LEO2_REED_SOL_HEADER_SHA256
             << "\",\"license\":\"BSD-3-Clause\",\"license_sha256\":\""
             << LEO2_JERASURE_LICENSE_SHA256
             << "\",\"gf_complete_commit\":\"" << LEO2_GF_COMPLETE_COMMIT
             << "\",\"gf_complete_tree\":\"" << LEO2_GF_COMPLETE_TREE
             << "\",\"gf_complete_header_sha256\":\""
             << LEO2_GF_COMPLETE_HEADER_SHA256
             << "\",\"gf_complete_license_sha256\":\""
             << LEO2_GF_COMPLETE_LICENSE_SHA256
             << "\",\"gf_complete_simd_flags\":\""
             << LEO2_GF_COMPLETE_SIMD_FLAGS
             << "\",\"field\":\"" << provider_field
             << "\",\"generator\":\"reed_sol_vandermonde_coding_matrix\""
             << ",\"wire_compatible\":false},\n"
             << "  \"parameters\":{\"K\":" << options.k << ",\"R\":" << options.r
             << ",\"requested_profile\":\"" << profile
             << "\",\"requested_field\":\"auto\",\"requested_backend\":\"auto\""
             << ",\"shard_bytes\":" << options.bytes << ",\"loss_count\":"
             << options.losses << ",\"missing_original_indices\":[";
        for (size_t i = 0; i < losses.size(); ++i)
        {
            if (i) json << ',';
            json << losses[i];
        }
        json << "],\"batch\":" << options.batch << ",\"reuse\":" << options.reuse
             << ",\"iterations\":" << options.iterations << ",\"warmup\":"
             << options.warmup << ",\"thread_count\":1,\"seed\":" << options.seed
             << "},\n  \"comparison_identity\":{\"leopard2_profile\":\"" << profile
             << "\",\"leopard2_parent\":" << parent << ",\"leopard2_field\":\""
             << leopard_field << "\",\"jerasure_field_advantage_from_padding\":"
             << ((width == 8 && parent > 256) ? "true" : "false")
             << ",\"scope\":\"public payload and repaired-output throughput only; field/basis representation, coordinates, generator matrices, and parity bytes differ\"},\n"
             << "  \"correctness\":{\"direct_source_round_trip\":true,"
             << "\"systematic_sources_immutable\":true,"
             << "\"independent_systematic_vandermonde_coefficients_checked\":"
             << independent.size() << ",\"independent_scalar_parity_mode\":\""
             << options.oracle << "\",\"independent_scalar_parity_checked_bytes_per_validation\":"
             << checked_parity << ",\"independent_scalar_parity_total_bytes_per_validation\":"
             << total_parity << ",\"independent_scalar_parity_validation_passes\":2},\n"
             << "  \"workload_digests\":{\"algorithm\":\"fnv1a64\","
             << "\"original_data\":\"" << HexU64(original_digest)
             << "\",\"transmitted_parity\":\"" << HexU64(parity_digest)
             << "\",\"recovered_originals\":\"" << HexU64(recovered_digest)
             << "\"},\n"
             << "  \"memory\":{\"alignment_bytes\":64,\"region_multiple_bytes\":8,"
             << "\"direct_application_buffers\":true,\"staging_copy_bytes_per_execution\":0,"
             << "\"encode_input_bytes_per_batch_call\":" << encode_input
             << ",\"encode_output_bytes_per_batch_call\":" << encode_output
             << ",\"decode_offered_bytes_per_batch_call\":" << offered
             << ",\"decode_selected_bytes_per_batch_call\":" << selected
             << ",\"decode_output_bytes_per_batch_call\":" << decode_output
             << ",\"encode_matrix_bytes\":" << encode.matrix.size() * sizeof(int)
             << ",\"decode_matrix_bytes\":" << decode.decoding_matrix.size() * sizeof(int)
             << ",\"decode_id_bytes\":" << decode.selected_ids.size() * sizeof(int)
             << "},\n  \"metrics\":{\"codec_setup\":";
        WriteSummary(json, codec_setup, 0, "", 0, "", false);
        json << ",\"encode_execution\":";
        WriteSummary(json, encode_execution, encode_input, "input_GB_per_s",
                     encode_output, "parity_output_GB_per_s", true);
        json << ",\"decode_plan_setup\":";
        WriteSummary(json, plan_setup, 0, "", 0, "", false);
        json << ",\"decode_execution\":";
        WriteSummary(json, decode_execution, offered, "offered_received_GB_per_s",
                     decode_output, "repaired_output_GB_per_s", true);
        const double amortized = decode_execution.median_us +
            plan_setup.median_us / static_cast<double>(options.reuse);
        json << ",\"decode_amortized_at_reuse\":{\"reuse_count\":" << options.reuse
             << ",\"derived_median_us_per_batch_call\":" << amortized
             << ",\"offered_received_GB_per_s\":";
        WriteRate(json, offered, amortized);
        json
             << ",\"repaired_output_GB_per_s\":";
        WriteRate(json, decode_output, amortized);
        json << "},\"rate_semantics\":\"offered_received counts all non-erased public shards for repair; Jerasure reads the recorded deterministic K-row subset; no-loss decode is a true no-op with null throughput\"}\n}\n";

        if (options.output == "-") std::cout << json.str();
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
    try { return Run(ParseOptions(argc, argv)); }
    catch (const std::exception& error)
    {
        std::cerr << "leopard2_jerasure_benchmark: " << error.what() << '\n';
        return 1;
    }
}
