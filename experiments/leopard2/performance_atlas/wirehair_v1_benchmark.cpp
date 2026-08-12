/*
    Copyright (c) 2017 Christopher A. Taylor.  All rights reserved.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the conditions in License.md
    are met.

    Standalone adapter for comparing the shipping Wirehair codec with the
    Leopard performance-atlas workload.  This file deliberately lives outside
    the default build: it is compiled against a pinned Wirehair checkout by
    the experiment runner and never adds a runtime dependency to Leopard.
*/

#include <wirehair/wirehair.h>
#include <gf256.h>

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
#include <new>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#if defined(_MSC_VER)
#include <malloc.h>
#endif

#ifndef LEO2_ATLAS_WIREHAIR_COMMIT
#error "LEO2_ATLAS_WIREHAIR_COMMIT must identify the Wirehair source"
#endif
#ifndef LEO2_ATLAS_PURE_AVX2
#error "LEO2_ATLAS_PURE_AVX2 must bind the adapter ISA profile"
#endif
#if LEO2_ATLAS_PURE_AVX2
#error "Shipping Wirehair contains target-qualified AVX-512 helpers"
#endif
#if !defined(LEO2_ATLAS_WIREHAIR_V1_COMPACT_PATH) || \
    !LEO2_ATLAS_WIREHAIR_V1_COMPACT_PATH
#error "Wirehair v1 atlas runs must force the compact XOR path"
#endif

namespace {

struct Options
{
    uint32_t k;
    uint32_t r;
    uint32_t block_bytes;
    uint32_t losses;
    size_t reuse;
    size_t iterations;
    size_t warmup;
    uint64_t seed;
    std::string output;

    Options()
        : k(224), r(32), block_bytes(1024), losses(1), reuse(8),
          iterations(7), warmup(2), seed(1), output("-")
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

class CodecOwner
{
public:
    CodecOwner() : codec_(NULL) {}
    ~CodecOwner() { wirehair_free(codec_); }

    WirehairCodec* Out()
    {
        wirehair_free(codec_);
        codec_ = NULL;
        return &codec_;
    }

    WirehairCodec Get() const { return codec_; }

private:
    CodecOwner(const CodecOwner&);
    CodecOwner& operator=(const CodecOwner&);
    WirehairCodec codec_;
};

class AlignedBytes
{
public:
    AlignedBytes() : data_(NULL), size_(0) {}
    ~AlignedBytes() { Release(); }

    void resize(size_t bytes)
    {
        Release();
        if (bytes == 0)
            return;
#if defined(_MSC_VER)
        data_ = static_cast<uint8_t*>(_aligned_malloc(bytes, 64));
#else
        void* pointer = NULL;
        if (posix_memalign(&pointer, 64, bytes) == 0)
            data_ = static_cast<uint8_t*>(pointer);
#endif
        if (!data_)
            throw std::bad_alloc();
        size_ = bytes;
        memset(data_, 0, bytes);
    }

    size_t size() const { return size_; }
    uint8_t& operator[](size_t index) { return data_[index]; }
    const uint8_t& operator[](size_t index) const { return data_[index]; }

private:
    AlignedBytes(const AlignedBytes&);
    AlignedBytes& operator=(const AlignedBytes&);

    void Release()
    {
#if defined(_MSC_VER)
        _aligned_free(data_);
#else
        free(data_);
#endif
        data_ = NULL;
        size_ = 0;
    }

    uint8_t* data_;
    size_t size_;
};

static void Fail(const std::string& message)
{
    throw std::runtime_error(message);
}

static void RequireCompactXorPath(const char* phase)
{
    // Target-qualified AVX-512 helpers in this Wirehair revision are reachable
    // only through its opt-in per-thread wide multi-source XOR mode.  The v1
    // codec normally leaves that mode disabled; make the requirement explicit
    // and fail if any benchmark phase returns with it enabled.
    if (gf256_set_thread_wide_xor(0) != 0)
        Fail(std::string("Wirehair wide XOR was active ") + phase);
}

static void Require(WirehairResult result, const char* operation)
{
    if (result != Wirehair_Success)
    {
        std::ostringstream message;
        message << operation << " failed: " << wirehair_result_string(result)
                << " (" << static_cast<int>(result) << ')';
        Fail(message.str());
    }
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

static uint32_t ParseUint32(const std::string& text, const char* option)
{
    const uint64_t value = ParseUnsigned(text, option);
    if (value > UINT32_MAX)
        Fail(std::string(option) + " exceeds uint32_t");
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
        Fail(std::string("missing value after ") + argv[index - 1]);
    return argv[index];
}

static void Usage(std::ostream& output, const char* program)
{
    output << "Usage: " << program << " [options]\n"
        << "  --k N             Source block count (2..64000)\n"
        << "  --r N             Generated repair blocks (default 32)\n"
        << "  --bytes N         Bytes per block\n"
        << "  --loss N          Random missing source blocks\n"
        << "  --reuse N         Steady encode calls per sample\n"
        << "  --iterations N    Timing samples\n"
        << "  --warmup N        Untimed repetitions\n"
        << "  --seed N          Deterministic data/pattern seed\n"
        << "  --json PATH       JSON output path, or - for stdout\n";
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
        else if (argument == "--k")
            options.k = ParseUint32(NeedValue(argc, argv, i), "--k");
        else if (argument == "--r")
            options.r = ParseUint32(NeedValue(argc, argv, i), "--r");
        else if (argument == "--bytes")
            options.block_bytes = ParseUint32(
                NeedValue(argc, argv, i), "--bytes");
        else if (argument == "--loss" || argument == "--losses")
            options.losses = ParseUint32(
                NeedValue(argc, argv, i), "--loss");
        else if (argument == "--reuse")
            options.reuse = ParseSize(
                NeedValue(argc, argv, i), "--reuse");
        else if (argument == "--iterations")
            options.iterations = ParseSize(
                NeedValue(argc, argv, i), "--iterations");
        else if (argument == "--warmup")
            options.warmup = ParseSize(
                NeedValue(argc, argv, i), "--warmup");
        else if (argument == "--seed")
            options.seed = ParseUnsigned(
                NeedValue(argc, argv, i), "--seed");
        else if (argument == "--json" || argument == "--output")
            options.output = NeedValue(argc, argv, i);
        else
            Fail("unknown argument: " + argument);
    }
    if (options.k < 2 || options.k > 64000)
        Fail("--k must be in 2..64000 for Wirehair v1");
    if (options.r == 0 || options.r > 1024)
        Fail("--r must be in 1..1024");
    if (options.block_bytes == 0 || options.block_bytes > INT32_MAX)
        Fail("--bytes must be in 1..INT32_MAX");
    if (options.losses == 0 || options.losses > options.k ||
        options.losses > options.r)
        Fail("--loss must be in 1..min(K,R)");
    if (options.reuse == 0 || options.iterations == 0)
        Fail("--reuse and --iterations must be positive");
    return options;
}

static size_t CheckedSize(uint64_t count, uint64_t bytes, const char* what)
{
    if (count != 0 && bytes > std::numeric_limits<size_t>::max() / count)
        Fail(std::string(what) + " size overflow");
    return static_cast<size_t>(count * bytes);
}

static std::vector<uint32_t> SelectLosses(const Options& options)
{
    std::vector<uint32_t> order(options.k);
    for (uint32_t i = 0; i < options.k; ++i)
        order[i] = i;
    XorShift64 random(options.seed ^ UINT64_C(0xd1b54a32d192ed03));
    for (size_t remaining = order.size(); remaining > 1; --remaining)
    {
        const size_t selected = static_cast<size_t>(
            random.Next() % remaining);
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
        deviations[i] = std::fabs(samples[i] - summary.median_us);
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

static double Rate(uint64_t bytes, double usec)
{
    return bytes == 0 || usec <= 0.0 ? 0.0 :
        static_cast<double>(bytes) / (usec * 1000.0);
}

static void WriteSummary(
    std::ostream& output,
    const Summary& summary,
    uint64_t input_bytes,
    const char* input_name,
    uint64_t output_bytes,
    const char* output_name,
    unsigned indent)
{
    const std::string spaces(indent, ' ');
    output << "{\n"
        << spaces << "  \"median_us\": " << summary.median_us << ",\n"
        << spaces << "  \"mad_us\": " << summary.mad_us << ",\n"
        << spaces << "  \"minimum_us\": " << summary.minimum_us << ",\n"
        << spaces << "  \"maximum_us\": " << summary.maximum_us << ",\n"
        << spaces << "  \"samples_us\": [";
    for (size_t i = 0; i < summary.samples_us.size(); ++i)
    {
        if (i != 0) output << ", ";
        output << summary.samples_us[i];
    }
    output << "],\n" << spaces << "  \"" << input_name << "\": "
        << Rate(input_bytes, summary.median_us) << ",\n"
        << spaces << "  \"" << output_name << "\": "
        << Rate(output_bytes, summary.median_us) << "\n"
        << spaces << '}';
}

static void WriteSetupSummary(
    std::ostream& output,
    const Summary& summary,
    unsigned indent)
{
    const std::string spaces(indent, ' ');
    output << "{\n"
        << spaces << "  \"median_us\": " << summary.median_us << ",\n"
        << spaces << "  \"mad_us\": " << summary.mad_us << ",\n"
        << spaces << "  \"minimum_us\": " << summary.minimum_us << ",\n"
        << spaces << "  \"maximum_us\": " << summary.maximum_us << ",\n"
        << spaces << "  \"samples_us\": [";
    for (size_t i = 0; i < summary.samples_us.size(); ++i)
    {
        if (i != 0) output << ", ";
        output << summary.samples_us[i];
    }
    output << "]\n" << spaces << '}';
}

struct Fixture
{
    AlignedBytes message;
    AlignedBytes repair;
    AlignedBytes recovered;
    std::vector<uint32_t> losses;
    uint32_t generated_repair_count;
    uint32_t consumed_repair_count;
    CodecOwner encoder;
};

static void GenerateRepairs(
    WirehairCodec encoder,
    const Options& options,
    uint32_t count,
    AlignedBytes& output)
{
    if (count > output.size() / options.block_bytes)
        Fail("repair output capacity is too small");
    for (uint32_t i = 0; i < count; ++i)
    {
        uint32_t bytes_out = 0;
        Require(wirehair_encode(encoder, options.k + i,
            &output[static_cast<size_t>(i) * options.block_bytes],
            options.block_bytes, &bytes_out), "wirehair_encode");
        if (bytes_out != options.block_bytes)
            Fail("Wirehair emitted a partial repair block");
    }
}

static uint32_t DecodeOnce(
    const Options& options,
    const Fixture& fixture,
    AlignedBytes& recovered,
    bool verify)
{
    CodecOwner decoder;
    Require(wirehair_decoder_create_ex(NULL,
        static_cast<uint64_t>(options.k) * options.block_bytes,
        options.block_bytes, decoder.Out()), "wirehair_decoder_create_ex");

    WirehairResult result = Wirehair_NeedMore;
    for (uint32_t i = 0; i < options.k && result == Wirehair_NeedMore; ++i)
    {
        if (IsLost(fixture.losses, i))
            continue;
        result = wirehair_decode(decoder.Get(), i,
            &fixture.message[static_cast<size_t>(i) * options.block_bytes],
            options.block_bytes);
        if (result != Wirehair_NeedMore && result != Wirehair_Success)
            Require(result, "wirehair_decode systematic");
    }

    uint32_t repairs = 0;
    while (result == Wirehair_NeedMore &&
           repairs < fixture.generated_repair_count)
    {
        result = wirehair_decode(decoder.Get(), options.k + repairs,
            &fixture.repair[
                static_cast<size_t>(repairs) * options.block_bytes],
            options.block_bytes);
        ++repairs;
        if (result != Wirehair_NeedMore && result != Wirehair_Success)
            Require(result, "wirehair_decode repair");
    }
    if (result != Wirehair_Success)
        Fail("Wirehair needed more than the retained repair budget");

    for (size_t i = 0; i < fixture.losses.size(); ++i)
    {
        uint32_t bytes_out = 0;
        uint8_t* destination = &recovered[i * options.block_bytes];
        Require(wirehair_recover_block_ex(decoder.Get(), fixture.losses[i],
            destination, options.block_bytes, &bytes_out),
            "wirehair_recover_block_ex");
        if (bytes_out != options.block_bytes)
            Fail("Wirehair recovered a partial source block");
        if (verify)
        {
            const uint8_t* expected = &fixture.message[
                static_cast<size_t>(fixture.losses[i]) * options.block_bytes];
            if (memcmp(expected, destination, options.block_bytes) != 0)
                Fail("Wirehair restored a source block incorrectly");
        }
    }
    return repairs;
}

static std::string Run(const Options& options)
{
    Require(wirehair_init(), "wirehair_init");
    gf256_x86_cpu_features active_features = {};
    gf256_get_active_x86_cpu_features(&active_features);
    if (!active_features.AVX2)
        Fail("Wirehair v1 AVX2 path requested but AVX2 is not active");
    if (active_features.GFNI)
        Fail("Wirehair GFNI is active despite the atlas compile flags");
    RequireCompactXorPath("before fixture construction");
    Fixture fixture;
    fixture.losses = SelectLosses(options);
    fixture.message.resize(CheckedSize(
        options.k, options.block_bytes, "message"));
    XorShift64 data_random(options.seed ^ UINT64_C(0x9e3779b97f4a7c15));
    for (size_t i = 0; i < fixture.message.size(); ++i)
        fixture.message[i] = static_cast<uint8_t>(data_random.Next() >> 56);

    fixture.generated_repair_count = options.r + 32;
    fixture.repair.resize(CheckedSize(fixture.generated_repair_count,
        options.block_bytes, "repair storage"));
    fixture.recovered.resize(CheckedSize(options.losses,
        options.block_bytes, "recovered storage"));
    Require(wirehair_encoder_create_ex(NULL, &fixture.message[0],
        fixture.message.size(), options.block_bytes, fixture.encoder.Out()),
        "wirehair_encoder_create_ex");
    GenerateRepairs(fixture.encoder.Get(), options,
        fixture.generated_repair_count, fixture.repair);
    fixture.consumed_repair_count = DecodeOnce(
        options, fixture, fixture.recovered, true);
    if (fixture.consumed_repair_count < options.losses)
        Fail("Wirehair consumed fewer repairs than missing source blocks");

    Fnv1a64 original_digest;
    Fnv1a64 parity_digest;
    Fnv1a64 recovered_digest;
    original_digest.Add(&fixture.message[0], fixture.message.size());
    parity_digest.Add(&fixture.repair[0], CheckedSize(
        options.r, options.block_bytes, "parity digest"));
    recovered_digest.Add(&fixture.recovered[0], fixture.recovered.size());

    AlignedBytes timed_repair;
    timed_repair.resize(CheckedSize(
        options.r, options.block_bytes, "timed repair"));
    AlignedBytes timed_recovered;
    timed_recovered.resize(fixture.recovered.size());
    const auto steady_encode = [&]() {
        GenerateRepairs(fixture.encoder.Get(), options, options.r, timed_repair);
    };
    const auto encoder_setup = [&]() {
        CodecOwner temporary;
        Require(wirehair_encoder_create_ex(NULL, &fixture.message[0],
            fixture.message.size(), options.block_bytes, temporary.Out()),
            "timed wirehair_encoder_create_ex");
    };
    const auto first_encode = [&]() {
        CodecOwner temporary;
        Require(wirehair_encoder_create_ex(NULL, &fixture.message[0],
            fixture.message.size(), options.block_bytes, temporary.Out()),
            "timed wirehair_encoder_create_ex");
        GenerateRepairs(temporary.Get(), options, options.r, timed_repair);
    };
    const auto first_decode = [&]() {
        const uint32_t used = DecodeOnce(
            options, fixture, timed_recovered, false);
        if (used != fixture.consumed_repair_count)
            Fail("Wirehair repair overhead changed between repetitions");
    };

    for (size_t i = 0; i < options.warmup; ++i)
    {
        steady_encode();
        first_encode();
        first_decode();
    }
    RequireCompactXorPath("after warmup");
    const Summary setup = Measure(options.iterations, 1, encoder_setup);
    RequireCompactXorPath("after setup timing");
    const Summary encode = Measure(
        options.iterations, options.reuse, steady_encode);
    RequireCompactXorPath("after encode timing");
    const Summary encode_first = Measure(
        options.iterations, options.reuse, first_encode);
    RequireCompactXorPath("after setup-inclusive encode timing");
    const Summary decode_first = Measure(
        options.iterations, options.reuse, first_decode);
    RequireCompactXorPath("after decode timing");

    if (memcmp(&timed_repair[0], &fixture.repair[0], timed_repair.size()) != 0)
        Fail("timed Wirehair repair output differs from correctness fixture");
    for (size_t i = 0; i < fixture.losses.size(); ++i)
    {
        const uint8_t* expected = &fixture.message[
            static_cast<size_t>(fixture.losses[i]) * options.block_bytes];
        const uint8_t* actual = &timed_recovered[i * options.block_bytes];
        if (memcmp(expected, actual, options.block_bytes) != 0)
            Fail("timed Wirehair recovery differs from source data");
    }

    const uint64_t source_bytes =
        static_cast<uint64_t>(options.k) * options.block_bytes;
    const uint64_t parity_bytes =
        static_cast<uint64_t>(options.r) * options.block_bytes;
    const uint64_t received_bytes =
        static_cast<uint64_t>(options.k - options.losses +
            fixture.consumed_repair_count) * options.block_bytes;
    const uint64_t repaired_bytes =
        static_cast<uint64_t>(options.losses) * options.block_bytes;

    std::ostringstream json;
    json << std::fixed << std::setprecision(6)
        << "{\n"
        << "  \"schema\": \"leopard2-performance-atlas-wirehair-v1/v2\",\n"
        << "  \"build\": {\"wirehair_source_commit\": \""
        << LEO2_ATLAS_WIREHAIR_COMMIT << "\", \"pure_avx2\": "
        << (LEO2_ATLAS_PURE_AVX2 ? "true" : "false")
        << ", \"wirehair_abi_version\": " << WIREHAIR_VERSION
        << ", \"wire_profile_version\": "
        << WIREHAIR_WIRE_PROFILE_VERSION
        << ", \"wire_profile_id\": " << WIREHAIR_LEGACY_PROFILE_CURRENT
        << ", \"isa_claim\": \"wirehair_v1_compact_path_avx2\""
        << ", \"target_qualified_avx512_helpers_present\": true"
        << ", \"wirehair_v1_wide_xor_forced_off\": true"
        << ", \"runtime_wide_xor_enabled\": false"
        << ", \"measured_path_avx512\": false"
        << ", \"active_x86_features\": {\"ssse3\": "
        << (active_features.SSSE3 ? "true" : "false")
        << ", \"avx2\": " << (active_features.AVX2 ? "true" : "false")
        << ", \"gfni\": " << (active_features.GFNI ? "true" : "false")
        << ", \"avx512\": "
        << (active_features.AVX512 ? "true" : "false") << "}"
        << ", \"cplusplus\": "
        << __cplusplus << "},\n"
        << "  \"parameters\": {\n"
        << "    \"K\": " << options.k << ", \"R\": " << options.r
        << ", \"shard_bytes\": " << options.block_bytes
        << ", \"loss_count\": " << options.losses << ",\n"
        << "    \"missing_original_indices\": ";
    WriteIndexArray(json, fixture.losses);
    json << ",\n    \"reuse\": " << options.reuse
        << ", \"iterations\": " << options.iterations
        << ", \"warmup\": " << options.warmup
        << ", \"seed\": " << options.seed
        << ", \"batch\": 1, \"thread_count\": 1\n"
        << "  },\n"
        << "  \"correctness\": {\"round_trip\": true},\n"
        << "  \"workload_digests\": {\n"
        << "    \"algorithm\": \"fnv1a64\",\n"
        << "    \"original_data\": \"" << original_digest.Hex() << "\",\n"
        << "    \"generated_repair\": \"" << parity_digest.Hex() << "\",\n"
        << "    \"recovered_originals\": \""
        << recovered_digest.Hex() << "\"\n"
        << "  },\n"
        << "  \"decode_input\": {\n"
        << "    \"surviving_systematic_blocks\": "
        << options.k - options.losses << ",\n"
        << "    \"repair_blocks_consumed\": "
        << fixture.consumed_repair_count << ",\n"
        << "    \"extra_repair_blocks\": "
        << fixture.consumed_repair_count - options.losses << ",\n"
        << "    \"arrival_order\": "
        << "\"surviving_systematic_ascending_then_repair_ascending\"\n"
        << "  },\n"
        << "  \"path_semantics\": {\n"
        << "    \"codec\": \"shipping_wirehair_v1\",\n"
        << "    \"threading\": \"single_caller_thread\",\n"
        << "    \"wide_xor\": \"forced_off_on_benchmark_thread\",\n"
        << "    \"avx512_target_helpers\": "
        << "\"present_but_unreachable_from_measured_v1_compact_path\"\n"
        << "  },\n"
        << "  \"timing_semantics\": {\n"
        << "    \"message_precode_setup\": "
        << "\"fresh wirehair_encoder_create_ex; no repair emission\",\n"
        << "    \"encode_execution\": "
        << "\"reuse one message-precode encoder; emit exactly R repairs; "
        << "exclude encoder creation\",\n"
        << "    \"encode_including_setup\": "
        << "\"fresh wirehair_encoder_create_ex then emit exactly R repairs\",\n"
        << "    \"decode_including_setup\": "
        << "\"fresh decoder; ingest surviving systematic blocks then "
        << "ascending repair blocks through solve; recover missing originals\"\n"
        << "  },\n"
        << "  \"metrics\": {\n"
        << "    \"message_precode_setup\": ";
    WriteSetupSummary(json, setup, 4);
    json << ",\n    \"encode_execution\": ";
    WriteSummary(json, encode, source_bytes, "message_equivalent_GB_per_s",
        parity_bytes, "repair_output_GB_per_s", 4);
    json << ",\n    \"encode_including_setup\": ";
    WriteSummary(json, encode_first, source_bytes,
        "message_equivalent_GB_per_s", parity_bytes,
        "repair_output_GB_per_s", 4);
    json << ",\n    \"decode_including_setup\": ";
    WriteSummary(json, decode_first, received_bytes,
        "received_input_GB_per_s", repaired_bytes,
        "repaired_output_GB_per_s", 4);
    json << "\n  }\n}\n";
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
        std::cerr << "Wirehair v1 atlas benchmark: " << error.what()
                  << std::endl;
        return 1;
    }
}
