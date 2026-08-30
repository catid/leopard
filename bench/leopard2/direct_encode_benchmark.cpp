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

#ifndef LEO2_ENABLE_TEST_HOOKS
#error "direct_encode_benchmark requires LEO2_ENABLE_TEST_HOOKS"
#endif
#if !defined(LEO2_BENCHMARK_SOURCE_ATTESTATION) || \
    !defined(LEO2_BENCHMARK_SOURCE_ATTESTATION_HEADER)
#error "direct_encode_benchmark requires generated source attestation"
#endif
#include LEO2_BENCHMARK_SOURCE_ATTESTATION_HEADER
#if !defined(LEO2_BENCHMARK_SOURCE_COMMIT) || \
    !defined(LEO2_BENCHMARK_SOURCE_TREE) || \
    !defined(LEO2_BENCHMARK_SOURCE_TRACKED_DIRTY) || \
    !defined(LEO2_BENCHMARK_BUILD_VARIANT) || \
    !defined(LEO2_BENCHMARK_BUILD_TYPE) || \
    !defined(LEO2_BENCHMARK_BUILD_CONFIGURATION_SHA256)
#error "direct benchmark requires source and build-variant attestation"
#endif

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
#include <locale>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#if defined(_MSC_VER)
#include <malloc.h>
#endif

namespace {

#define LEO2_STRINGIZE_DETAIL(value) #value
#define LEO2_STRINGIZE(value) LEO2_STRINGIZE_DETAIL(value)

enum EncodeMode
{
    MODE_UNSPECIFIED,
    MODE_AUTO,
    MODE_DIRECT,
    MODE_TRANSFORM
};

struct Options
{
    uint32_t k;
    uint32_t r;
    leo2_profile profile;
    leo2_field field;
    uint64_t bytes;
    uint32_t q;
    bool q_was_set;
    std::string requested_parity;
    bool requested_parity_was_set;
    size_t batch;
    size_t reuse;
    size_t iterations;
    size_t warmups;
    uint32_t threads;
    uint64_t seed;
    EncodeMode mode;
    std::string output;

    Options()
        : k(8)
        , r(8)
        , profile(LEO2_PROFILE_AUTO)
        , field(LEO2_FIELD_AUTO)
        , bytes(1024)
        , q(0)
        , q_was_set(false)
        , requested_parity("all")
        , requested_parity_was_set(false)
        , batch(1)
        , reuse(64)
        , iterations(15)
        , warmups(4)
        , threads(1)
        , seed(1)
        , mode(MODE_UNSPECIFIED)
        , output("-")
    {}
};

struct Summary
{
    double median_us;
    double mad_us;
    double minimum_us;
    double maximum_us;
};

static void Fail(const std::string& message)
{
    throw std::runtime_error(message);
}

static void RequireLeo2(leo2_result result, const char* operation)
{
    if (result != LEO2_SUCCESS)
    {
        std::ostringstream stream;
        stream << operation << " failed: " << leo2_result_string(result)
               << " (" << static_cast<int>(result) << ')';
        Fail(stream.str());
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
    if (value > std::numeric_limits<uint32_t>::max())
        Fail(std::string("value for ") + option + " exceeds uint32");
    return static_cast<uint32_t>(value);
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
    size_t digit_count = 0;
    while (digit_count < text.size() && text[digit_count] >= '0' &&
           text[digit_count] <= '9')
        ++digit_count;
    if (digit_count == 0)
        Fail("invalid value for --bytes: " + text);

    const uint64_t base = ParseUnsigned(text.substr(0, digit_count), "--bytes");
    std::string suffix = text.substr(digit_count);
    std::transform(suffix.begin(), suffix.end(), suffix.begin(),
        [](unsigned char character) {
            return static_cast<char>(std::tolower(character));
        });

    uint64_t multiplier = 1;
    if (suffix.empty() || suffix == "b")
        multiplier = 1;
    else if (suffix == "k" || suffix == "kb")
        multiplier = UINT64_C(1000);
    else if (suffix == "m" || suffix == "mb")
        multiplier = UINT64_C(1000000);
    else if (suffix == "g" || suffix == "gb")
        multiplier = UINT64_C(1000000000);
    else if (suffix == "ki" || suffix == "kib")
        multiplier = UINT64_C(1024);
    else if (suffix == "mi" || suffix == "mib")
        multiplier = UINT64_C(1024) * UINT64_C(1024);
    else if (suffix == "gi" || suffix == "gib")
        multiplier = UINT64_C(1024) * UINT64_C(1024) * UINT64_C(1024);
    else
        Fail("invalid suffix for --bytes: " + suffix);
    if (base > std::numeric_limits<uint64_t>::max() / multiplier)
        Fail("--bytes overflows uint64");
    return base * multiplier;
}

static leo2_profile ParseProfile(const std::string& text)
{
    if (text == "auto") return LEO2_PROFILE_AUTO;
    if (text == "high" || text == "legacy-high" ||
        text == "legacy-high-v1" || text == "legacy_high_v1")
        return LEO2_PROFILE_LEGACY_HIGH_V1;
    if (text == "low" || text == "low-v1" || text == "low_v1")
        return LEO2_PROFILE_LOW_V1;
    Fail("invalid --profile: " + text);
    return LEO2_PROFILE_AUTO;
}

static leo2_field ParseField(const std::string& text)
{
    if (text == "auto") return LEO2_FIELD_AUTO;
    if (text == "gf8") return LEO2_FIELD_GF8;
    if (text == "gf16") return LEO2_FIELD_GF16;
    Fail("invalid --field: " + text);
    return LEO2_FIELD_AUTO;
}

static const char* NeedValue(int argc, char** argv, int& index)
{
    if (++index >= argc)
        Fail(std::string("missing value after ") + argv[index - 1]);
    return argv[index];
}

static void SelectMode(Options& options, EncodeMode mode, const char* option)
{
    if (options.mode != MODE_UNSPECIFIED && options.mode != mode)
        Fail(std::string(option) + " conflicts with the selected encode mode");
    options.mode = mode;
}

static void Usage(std::ostream& output, const char* program)
{
    output
        << "Usage: " << program
        << " --auto|--force-direct|--force-transform [options]\n"
        << "  --k N                       Original shard count (default 8)\n"
        << "  --r N                       Recovery shard count (default 8)\n"
        << "  --profile NAME              auto, high, low (default auto)\n"
        << "  --field NAME                auto, gf8, gf16 (default auto)\n"
        << "  --bytes N[KiB|MiB]          Bytes per shard (default 1KiB)\n"
        << "  --q N                       Requested parity count; without a mask,\n"
        << "                              request prefix [0,N) (default R)\n"
        << "  --requested-parity MASK     all, or comma-separated indices/ranges,\n"
        << "                              for example 0,2-3,7\n"
        << "  --batch N                   Independent stripes per call (default 1)\n"
        << "  --reuse N                   Calls per timing sample (default 64)\n"
        << "  --iterations N              Timing samples (default 15)\n"
        << "  --warmups N                 Untimed calls/setups (default 4)\n"
        << "  --threads N                 Context thread count (default 1)\n"
        << "  --seed N                    Deterministic data seed (default 1)\n"
        << "  --auto                      Time ordinary production AUTO dispatch\n"
        << "  --force-direct              Time the bounded direct encoder\n"
        << "  --force-transform           Time the transform encoder\n"
        << "  --mode auto|direct|transform  Equivalent explicit mode spelling\n"
        << "  --json PATH                 JSON output path, or - (default -)\n"
        << "  --help                      Show this message\n";
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
        else if (argument == "--profile")
            options.profile = ParseProfile(NeedValue(argc, argv, i));
        else if (argument == "--field")
            options.field = ParseField(NeedValue(argc, argv, i));
        else if (argument == "--bytes")
            options.bytes = ParseBytes(NeedValue(argc, argv, i));
        else if (argument == "--q")
        {
            options.q = ParseUint32(NeedValue(argc, argv, i), "--q");
            options.q_was_set = true;
        }
        else if (argument == "--requested-parity" ||
                 argument == "--requested-parity-mask" ||
                 argument == "--requested-mask")
        {
            options.requested_parity = NeedValue(argc, argv, i);
            options.requested_parity_was_set = true;
        }
        else if (argument == "--batch")
            options.batch = ParseSize(NeedValue(argc, argv, i), "--batch");
        else if (argument == "--reuse")
            options.reuse = ParseSize(NeedValue(argc, argv, i), "--reuse");
        else if (argument == "--iterations")
            options.iterations = ParseSize(
                NeedValue(argc, argv, i), "--iterations");
        else if (argument == "--warmup" || argument == "--warmups")
            options.warmups = ParseSize(NeedValue(argc, argv, i), "--warmups");
        else if (argument == "--threads" || argument == "--thread-count")
            options.threads = ParseUint32(NeedValue(argc, argv, i), "--threads");
        else if (argument == "--seed")
            options.seed = ParseUnsigned(NeedValue(argc, argv, i), "--seed");
        else if (argument == "--auto")
            SelectMode(options, MODE_AUTO, "--auto");
        else if (argument == "--force-direct")
            SelectMode(options, MODE_DIRECT, "--force-direct");
        else if (argument == "--force-transform")
            SelectMode(options, MODE_TRANSFORM, "--force-transform");
        else if (argument == "--mode")
        {
            const std::string mode = NeedValue(argc, argv, i);
            if (mode == "auto")
                SelectMode(options, MODE_AUTO, "--mode auto");
            else if (mode == "direct")
                SelectMode(options, MODE_DIRECT, "--mode direct");
            else if (mode == "transform")
                SelectMode(options, MODE_TRANSFORM, "--mode transform");
            else
                Fail("invalid --mode: " + mode);
        }
        else if (argument == "--json" || argument == "--output")
            options.output = NeedValue(argc, argv, i);
        else
            Fail("unknown argument: " + argument);
    }

    if (options.k == 0 || options.r == 0)
        Fail("--k and --r must be positive");
    if (options.bytes == 0 || options.bytes > std::numeric_limits<size_t>::max())
        Fail("--bytes must be positive and fit in size_t");
    if (options.batch == 0 || options.reuse == 0 || options.iterations == 0 ||
        options.threads == 0)
        Fail("--batch, --reuse, --iterations, and --threads must be positive");
    if (options.mode == MODE_UNSPECIFIED)
        Fail("select exactly one encode mode");
    return options;
}

static std::vector<uint8_t> ParseRequestedMask(Options& options)
{
    std::vector<uint8_t> requested(options.r, 0);
    if (!options.requested_parity_was_set)
    {
        const uint32_t q = options.q_was_set ? options.q : options.r;
        if (q == 0 || q > options.r)
            Fail("--q must be in [1,R]");
        for (uint32_t i = 0; i < q; ++i)
            requested[i] = 1;
        options.q = q;
        return requested;
    }

    const std::string& mask = options.requested_parity;
    if (mask == "all")
    {
        std::fill(requested.begin(), requested.end(), 1);
    }
    else if (mask.empty() || mask == "none")
    {
        Fail("--requested-parity must select at least one parity shard");
    }
    else
    {
        size_t begin = 0;
        while (begin < mask.size())
        {
            const size_t comma = mask.find(',', begin);
            const size_t end = comma == std::string::npos ? mask.size() : comma;
            if (end == begin)
                Fail("empty item in --requested-parity: " + mask);
            const std::string item = mask.substr(begin, end - begin);
            const size_t dash = item.find('-');
            uint32_t first = 0;
            uint32_t last = 0;
            if (dash == std::string::npos)
            {
                first = last = ParseUint32(item, "--requested-parity");
            }
            else
            {
                if (dash == 0 || dash + 1 == item.size() ||
                    item.find('-', dash + 1) != std::string::npos)
                    Fail("invalid range in --requested-parity: " + item);
                first = ParseUint32(
                    item.substr(0, dash), "--requested-parity");
                last = ParseUint32(
                    item.substr(dash + 1), "--requested-parity");
                if (last < first)
                    Fail("descending range in --requested-parity: " + item);
            }
            if (last >= options.r)
                Fail("--requested-parity index is outside [0,R)");
            for (uint32_t index = first; index <= last; ++index)
            {
                if (requested[index])
                    Fail("duplicate index in --requested-parity");
                requested[index] = 1;
            }
            if (comma == std::string::npos)
                break;
            begin = comma + 1;
        }
    }

    const uint32_t q = static_cast<uint32_t>(
        std::count(requested.begin(), requested.end(), static_cast<uint8_t>(1)));
    if (q == 0)
        Fail("--requested-parity must select at least one parity shard");
    if (options.q_was_set && options.q != q)
        Fail("--q does not equal the requested parity mask cardinality");
    options.q = q;
    return requested;
}

static void AlignedFree(void* pointer)
{
#if defined(_MSC_VER)
    _aligned_free(pointer);
#else
    free(pointer);
#endif
}

static void* AlignedAllocate(size_t alignment, size_t bytes)
{
#if defined(_MSC_VER)
    return _aligned_malloc(bytes, alignment);
#else
    void* pointer = NULL;
    return posix_memalign(&pointer, alignment, bytes) == 0 ? pointer : NULL;
#endif
}

class AlignedBuffer
{
public:
    AlignedBuffer()
        : data_(NULL)
        , size_(0)
    {}

    ~AlignedBuffer()
    {
        AlignedFree(data_);
    }

    void Reset(size_t bytes)
    {
        AlignedFree(data_);
        data_ = NULL;
        size_ = 0;
        if (bytes == 0)
            return;
        const size_t alignment = std::max<size_t>(64, leo2_scratch_alignment());
        data_ = AlignedAllocate(alignment, bytes);
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

static size_t CheckedSize(uint64_t count, uint64_t bytes, const char* what)
{
    if (count != 0 && bytes > std::numeric_limits<size_t>::max() / count)
        Fail(std::string(what) + " size overflows size_t");
    return static_cast<size_t>(count * bytes);
}

static uint64_t CheckedProduct(uint64_t left, uint64_t right, const char* what)
{
    if (left != 0 && right > std::numeric_limits<uint64_t>::max() / left)
        Fail(std::string(what) + " overflows uint64");
    return left * right;
}

struct Stripe
{
    AlignedBuffer original_storage;
    AlignedBuffer direct_storage;
    AlignedBuffer transform_storage;
    AlignedBuffer scratch;
    std::vector<const void*> original;
    std::vector<void*> direct_recovery;
    std::vector<void*> transform_recovery;
};

static void InitializeStripe(
    Stripe& stripe,
    const Options& options,
    const std::vector<uint8_t>& requested,
    size_t scratch_bytes,
    size_t stripe_index)
{
    const size_t original_bytes = CheckedSize(
        options.k, options.bytes, "original storage");
    const size_t recovery_bytes = CheckedSize(
        options.r, options.bytes, "recovery storage");
    stripe.original_storage.Reset(original_bytes);
    stripe.direct_storage.Reset(recovery_bytes);
    stripe.transform_storage.Reset(recovery_bytes);
    stripe.scratch.Reset(scratch_bytes);
    stripe.original.resize(options.k);
    stripe.direct_recovery.assign(options.r, static_cast<void*>(NULL));
    stripe.transform_recovery.assign(options.r, static_cast<void*>(NULL));

    XorShift64 random(options.seed ^
        (UINT64_C(0x9e3779b97f4a7c15) * (stripe_index + 1)));
    for (size_t i = 0; i < original_bytes; ++i)
        stripe.original_storage.bytes()[i] =
            static_cast<uint8_t>(random.Next() >> 56);
    memset(stripe.direct_storage.data(), 0xa5, recovery_bytes);
    memset(stripe.transform_storage.data(), 0xa5, recovery_bytes);

    for (uint32_t i = 0; i < options.k; ++i)
        stripe.original[i] = stripe.original_storage.bytes() +
            static_cast<size_t>(i) * static_cast<size_t>(options.bytes);
    for (uint32_t i = 0; i < options.r; ++i)
    {
        if (!requested[i])
            continue;
        const size_t offset =
            static_cast<size_t>(i) * static_cast<size_t>(options.bytes);
        stripe.direct_recovery[i] = stripe.direct_storage.bytes() + offset;
        stripe.transform_recovery[i] = stripe.transform_storage.bytes() + offset;
    }
}

static leo2_test_encode_mode TestMode(EncodeMode mode)
{
    if (mode == MODE_AUTO)
        return LEO2_TEST_ENCODE_AUTO;
    if (mode == MODE_DIRECT)
        return LEO2_TEST_ENCODE_FORCE_DIRECT;
    return LEO2_TEST_ENCODE_FORCE_TRANSFORM;
}

static const char* ModeName(EncodeMode mode)
{
    if (mode == MODE_AUTO)
        return "auto";
    return mode == MODE_DIRECT ? "force_direct" : "force_transform";
}

#if defined(LEO2_EXPERIMENT_HIGH_DIRECT_ENCODE) && \
    LEO2_EXPERIMENT_HIGH_DIRECT_ENCODE_AUTO
static bool ExpectedMeasuredHighFullDirect(const Options& options)
{
    const bool exact_k5_r5 = options.k == 5 && options.r == 5;
    if (exact_k5_r5)
    {
        if (options.bytes <= 14)
            return true;
        if (options.bytes >= 16 && options.bytes <= 60 &&
            options.bytes % 4 != 3)
            return true;
        if (options.bytes >= 64 && options.bytes <= 70)
            return true;
        if (options.bytes == 96)
            return true;
        return options.bytes >= 128 && options.bytes <= 4096 &&
            (options.bytes & 63U) == 0;
    }

    switch (options.bytes)
    {
    case 1: case 2:
    case 4: case 5: case 6:
    case 8: case 9: case 10:
    case 12: case 13: case 14:
    case 16: case 17: case 18:
    case 20: case 21: case 22:
    case 24: case 25: case 26:
    case 28:
    case 32: case 33: case 34:
    case 36: case 38:
    case 40: case 41:
    case 48:
    case 65:
        return true;
    case 3:
    case 64:
        return options.k != 16 || options.r != 8;
    default:
        return false;
    }
}
#endif

static bool ExpectedAutoDirectPath(
    const Options& options,
    const leo2_codec* codec,
    const leo2_context* context,
    bool direct_capable)
{
    const leo2_profile profile = leo2_codec_profile(codec);
    if (!direct_capable || options.k < 2)
        return false;
    const leo2_backend backend = leo2_context_backend(context);
    if (profile == LEO2_PROFILE_LEGACY_HIGH_V1)
    {
        if (leo2_codec_field(codec) != LEO2_FIELD_GF8 ||
            backend != LEO2_BACKEND_AVX2)
            return false;
#if defined(LEO2_EXPERIMENT_HIGH_DIRECT_ENCODE) && \
    LEO2_EXPERIMENT_HIGH_DIRECT_ENCODE_AUTO
        if (options.k >= 5 && options.k <= 16 &&
            options.r >= 5 && options.r <= 8 &&
            options.q == options.r &&
            ExpectedMeasuredHighFullDirect(options))
            return true;
#endif
#if (defined(LEO2_EXPERIMENT_HIGH_DIRECT_ENCODE) && \
     LEO2_EXPERIMENT_HIGH_DIRECT_ENCODE_AUTO) || \
    defined(LEO2_EXPERIMENT_HIGH_SPARSE_DIRECT_ENCODE)
        return options.q == 1 && options.bytes >= 1024 &&
            (options.bytes & 63U) == 0 && options.r > 1;
#else
        return false;
#endif
    }
    if (options.q != 1 || options.bytes < 1024 ||
        (options.bytes & 63U) != 0)
        return false;
    if (profile != LEO2_PROFILE_LOW_V1)
        return false;
    if (backend == LEO2_BACKEND_SCALAR)
        return options.k >= 3;
    return backend == LEO2_BACKEND_SSSE3 ||
        backend == LEO2_BACKEND_AVX2 ||
        backend == LEO2_BACKEND_AVX512 ||
        backend == LEO2_BACKEND_GFNI;
}

static const char* ProfileName(leo2_profile profile)
{
    switch (profile)
    {
    case LEO2_PROFILE_AUTO: return "auto";
    case LEO2_PROFILE_LEGACY_HIGH_V1: return "legacy_high_v1";
    case LEO2_PROFILE_LOW_V1: return "low_v1";
    case LEO2_PROFILE_EXACT_EXPERIMENTAL_V1: return "exact_experimental_v1";
    }
    return "unknown";
}

static const char* FieldName(leo2_field field)
{
    switch (field)
    {
    case LEO2_FIELD_AUTO: return "auto";
    case LEO2_FIELD_GF8: return "gf8";
    case LEO2_FIELD_GF16: return "gf16";
    }
    return "unknown";
}

static const char* BackendName(leo2_backend backend)
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

static const char* CompilerName()
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
    return LEO2_STRINGIZE(_MSC_FULL_VER);
#else
    return "unknown";
#endif
}

static leo2_codec* CreateCodec(
    leo2_context* context,
    const Options& options,
    EncodeMode mode)
{
    leo2_codec_options codec_options = {};
    codec_options.struct_size = sizeof(codec_options);
    leo2_codec* codec = NULL;
    RequireLeo2(leo2_codec_create(context, options.k, options.r,
        options.profile, options.field, &codec_options, &codec),
        "codec create");
    const leo2_result mode_result = leo2_test_codec_set_encode_mode(
        codec, TestMode(mode));
    if (mode_result != LEO2_SUCCESS)
    {
        leo2_codec_destroy(codec);
        RequireLeo2(mode_result, "set forced encode mode");
    }
    return codec;
}

static double Median(std::vector<double> values)
{
    const size_t middle = values.size() / 2;
    std::nth_element(values.begin(), values.begin() + middle, values.end());
    const double upper = values[middle];
    if ((values.size() & 1u) != 0)
        return upper;
    std::nth_element(
        values.begin(), values.begin() + middle - 1, values.begin() + middle);
    return (values[middle - 1] + upper) * 0.5;
}

static Summary Summarize(const std::vector<double>& samples)
{
    Summary summary;
    summary.median_us = Median(samples);
    std::vector<double> deviations(samples.size());
    for (size_t i = 0; i < samples.size(); ++i)
        deviations[i] = std::fabs(samples[i] - summary.median_us);
    summary.mad_us = Median(deviations);
    summary.minimum_us = *std::min_element(samples.begin(), samples.end());
    summary.maximum_us = *std::max_element(samples.begin(), samples.end());
    return summary;
}

static Summary MeasureCodecSetup(
    leo2_context* context,
    const Options& options)
{
    typedef std::chrono::steady_clock Clock;
    for (size_t i = 0; i < options.warmups; ++i)
    {
        leo2_codec* codec = CreateCodec(context, options, options.mode);
        leo2_codec_destroy(codec);
    }

    std::vector<double> samples;
    samples.reserve(options.iterations);
    for (size_t i = 0; i < options.iterations; ++i)
    {
        const Clock::time_point begin = Clock::now();
        leo2_codec* codec = CreateCodec(context, options, options.mode);
        leo2_codec_destroy(codec);
        const Clock::time_point end = Clock::now();
        samples.push_back(std::chrono::duration_cast<
            std::chrono::duration<double, std::micro> >(end - begin).count());
    }
    return Summarize(samples);
}

static Summary MeasureEncode(
    const leo2_codec* codec,
    const std::vector<leo2_encode_batch_item>& items,
    const Options& options)
{
    typedef std::chrono::steady_clock Clock;
    for (size_t i = 0; i < options.warmups; ++i)
        RequireLeo2(leo2_encode_batch(codec, &items[0], items.size()),
            "encode warmup");

    // samples is reserved before timing.  The batch items, shard storage, and
    // per-stripe scratch are likewise immutable/preallocated for every timed
    // call; errors are the only path that constructs an exception message.
    std::vector<double> samples;
    samples.reserve(options.iterations);
    for (size_t sample = 0; sample < options.iterations; ++sample)
    {
        const Clock::time_point begin = Clock::now();
        for (size_t call = 0; call < options.reuse; ++call)
            RequireLeo2(leo2_encode_batch(codec, &items[0], items.size()),
                "timed encode");
        const Clock::time_point end = Clock::now();
        const double microseconds = std::chrono::duration_cast<
            std::chrono::duration<double, std::micro> >(end - begin).count();
        samples.push_back(microseconds / static_cast<double>(options.reuse));
    }
    return Summarize(samples);
}

static void VerifyUntouched(
    const uint8_t* shard,
    size_t bytes,
    const char* path,
    size_t stripe,
    uint32_t recovery)
{
    for (size_t i = 0; i < bytes; ++i)
    {
        if (shard[i] != 0xa5)
        {
            std::ostringstream stream;
            stream << path << " modified unrequested parity " << recovery
                   << " in stripe " << stripe << " at byte " << i;
            Fail(stream.str());
        }
    }
}

static void VerifyParity(
    const std::vector<std::unique_ptr<Stripe> >& stripes,
    const Options& options,
    const std::vector<uint8_t>& requested)
{
    const size_t bytes = static_cast<size_t>(options.bytes);
    for (size_t stripe_index = 0; stripe_index < stripes.size(); ++stripe_index)
    {
        const Stripe& stripe = *stripes[stripe_index];
        for (uint32_t recovery = 0; recovery < options.r; ++recovery)
        {
            const size_t offset = static_cast<size_t>(recovery) * bytes;
            const uint8_t* direct = stripe.direct_storage.bytes() + offset;
            const uint8_t* transform = stripe.transform_storage.bytes() + offset;
            if (requested[recovery])
            {
                if (memcmp(direct, transform, bytes) != 0)
                {
                    std::ostringstream stream;
                    stream << "selected and transform-reference parity differ at stripe "
                           << stripe_index << ", recovery " << recovery;
                    Fail(stream.str());
                }
            }
            else
            {
                VerifyUntouched(direct, bytes, "direct", stripe_index, recovery);
                VerifyUntouched(
                    transform, bytes, "transform", stripe_index, recovery);
            }
        }
    }
}

static uint64_t ParityChecksum(
    const std::vector<std::unique_ptr<Stripe> >& stripes,
    const Options& options,
    const std::vector<uint8_t>& requested)
{
    uint64_t hash = UINT64_C(14695981039346656037);
    const size_t bytes = static_cast<size_t>(options.bytes);
    for (size_t stripe_index = 0; stripe_index < stripes.size(); ++stripe_index)
    {
        const Stripe& stripe = *stripes[stripe_index];
        for (uint32_t recovery = 0; recovery < options.r; ++recovery)
        {
            if (!requested[recovery])
                continue;
            const uint8_t* shard = stripe.transform_storage.bytes() +
                static_cast<size_t>(recovery) * bytes;
            for (size_t i = 0; i < bytes; ++i)
            {
                hash ^= shard[i];
                hash *= UINT64_C(1099511628211);
            }
        }
    }
    return hash;
}

static std::string JsonEscape(const std::string& input)
{
    std::ostringstream output;
    output.imbue(std::locale::classic());
    for (size_t i = 0; i < input.size(); ++i)
    {
        const unsigned char character = static_cast<unsigned char>(input[i]);
        switch (character)
        {
        case '\\': output << "\\\\"; break;
        case '"': output << "\\\""; break;
        case '\b': output << "\\b"; break;
        case '\f': output << "\\f"; break;
        case '\n': output << "\\n"; break;
        case '\r': output << "\\r"; break;
        case '\t': output << "\\t"; break;
        default:
            if (character < 0x20)
                output << "\\u" << std::hex << std::setw(4)
                       << std::setfill('0') << static_cast<unsigned>(character)
                       << std::dec;
            else
                output << input[i];
        }
    }
    return output.str();
}

static double GigabytesPerSecond(uint64_t bytes, double microseconds)
{
    if (bytes == 0 || microseconds <= 0)
        return 0;
    return static_cast<double>(bytes) / (microseconds * 1000.0);
}

static void WriteSummary(
    std::ostream& output,
    const Summary& summary,
    const char* suffix,
    unsigned indent)
{
    const std::string spaces(indent, ' ');
    output << "{\n"
           << spaces << "  \"median_us" << suffix << "\": "
           << summary.median_us << ",\n"
           << spaces << "  \"mad_us" << suffix << "\": "
           << summary.mad_us << ",\n"
           << spaces << "  \"minimum_us" << suffix << "\": "
           << summary.minimum_us << ",\n"
           << spaces << "  \"maximum_us" << suffix << "\": "
           << summary.maximum_us << "\n"
           << spaces << '}';
}

static void WriteEncodeSummary(
    std::ostream& output,
    const Summary& summary,
    uint64_t input_bytes,
    uint64_t output_bytes,
    unsigned indent)
{
    const std::string spaces(indent, ' ');
    output << "{\n"
           << spaces << "  \"median_us_per_batch_call\": "
           << summary.median_us << ",\n"
           << spaces << "  \"mad_us_per_batch_call\": "
           << summary.mad_us << ",\n"
           << spaces << "  \"minimum_us_per_batch_call\": "
           << summary.minimum_us << ",\n"
           << spaces << "  \"maximum_us_per_batch_call\": "
           << summary.maximum_us << ",\n"
           << spaces << "  \"logical_input_GB_per_s\": "
           << GigabytesPerSecond(input_bytes, summary.median_us) << ",\n"
           << spaces << "  \"requested_parity_output_GB_per_s\": "
           << GigabytesPerSecond(output_bytes, summary.median_us) << ",\n"
           << spaces << "  \"logical_input_plus_output_GB_per_s\": "
           << GigabytesPerSecond(
                CheckedProduct(1, input_bytes + output_bytes, "logical I/O"),
                summary.median_us) << "\n"
           << spaces << '}';
}

static int Run(Options options)
{
    leo2_context_options context_options = {};
    context_options.struct_size = sizeof(context_options);
    context_options.backend = LEO2_BACKEND_AUTO;
    context_options.thread_count = options.threads;
    leo2_context* context = NULL;
    RequireLeo2(leo2_context_create(&context_options, &context),
        "context create");

    leo2_codec* codec = NULL;
    try
    {
        codec = CreateCodec(context, options, options.mode);
        const bool direct_capable =
            leo2_test_codec_direct_encode_capable(codec) != 0;
        if (options.mode == MODE_DIRECT && !direct_capable)
            Fail("the exact cell is outside the bounded direct encoder domain");
        const std::vector<uint8_t> requested = ParseRequestedMask(options);

        size_t scratch_bytes = 0;
        RequireLeo2(leo2_encode_scratch_size(
            codec, options.bytes, &scratch_bytes), "encode scratch query");
        RequireLeo2(leo2_test_codec_set_encode_mode(
            codec, LEO2_TEST_ENCODE_FORCE_TRANSFORM),
            "select reference transform for scratch query");
        size_t reference_scratch_bytes = 0;
        RequireLeo2(leo2_encode_scratch_size(codec, options.bytes,
            &reference_scratch_bytes), "reference transform scratch query");
        RequireLeo2(leo2_test_codec_set_encode_mode(codec,
            TestMode(options.mode)), "restore selected encode mode");
        const size_t allocated_scratch_bytes =
            std::max(scratch_bytes, reference_scratch_bytes);

        std::vector<std::unique_ptr<Stripe> > stripes;
        stripes.reserve(options.batch);
        std::vector<leo2_encode_batch_item> direct_items(options.batch);
        std::vector<leo2_encode_batch_item> transform_items(options.batch);
        for (size_t i = 0; i < options.batch; ++i)
        {
            stripes.push_back(std::unique_ptr<Stripe>(new Stripe));
            InitializeStripe(
                *stripes.back(), options, requested,
                allocated_scratch_bytes, i);

            leo2_encode_batch_item& direct = direct_items[i];
            direct.shard_bytes = options.bytes;
            direct.original = &stripes.back()->original[0];
            direct.recovery = &stripes.back()->direct_recovery[0];
            direct.scratch = stripes.back()->scratch.data();
            direct.scratch_bytes = stripes.back()->scratch.size();

            leo2_encode_batch_item& transform = transform_items[i];
            transform.shard_bytes = options.bytes;
            transform.original = &stripes.back()->original[0];
            transform.recovery = &stripes.back()->transform_recovery[0];
            transform.scratch = stripes.back()->scratch.data();
            transform.scratch_bytes = stripes.back()->scratch.size();
        }

        // Keep a mature transform reference in direct_storage.  The selected
        // path writes transform_storage, so AUTO can be timed for arbitrary
        // K/R while retaining parity and unrequested-output checks.  For a
        // forced-direct cell this is the same direct-versus-transform oracle
        // as before, with the two storage labels merely reversed.
        int selected_direct = -1;
        RequireLeo2(leo2_test_codec_set_encode_mode(
            codec, LEO2_TEST_ENCODE_FORCE_TRANSFORM),
            "force transform validation");
        RequireLeo2(leo2_test_codec_encode_path(
            codec, options.bytes, options.q, &selected_direct),
            "query forced transform path");
        if (selected_direct != 0)
            Fail("forced transform path introspection selected direct");
        RequireLeo2(leo2_encode_batch(
            codec, &direct_items[0], direct_items.size()),
            "reference transform correctness encode");

        RequireLeo2(leo2_test_codec_set_encode_mode(codec,
            TestMode(options.mode)), "select timed encode path");
        RequireLeo2(leo2_test_codec_encode_path(
            codec, options.bytes, options.q, &selected_direct),
            "query timed encode path");
        const int expected_direct = options.mode == MODE_DIRECT ? 1 :
            (options.mode == MODE_TRANSFORM ? 0 :
             (ExpectedAutoDirectPath(
                  options, codec, context, direct_capable) ? 1 : 0));
        if (selected_direct != expected_direct)
            Fail("timed encode path introspection disagrees with selected mode");
        RequireLeo2(leo2_encode_batch(
            codec, &transform_items[0], transform_items.size()),
            "selected-path correctness encode");
        VerifyParity(stripes, options, requested);

        const Summary setup = MeasureCodecSetup(context, options);
        const Summary encode = MeasureEncode(codec, transform_items, options);

        // Re-check after timed execution so the reported checksum and
        // correctness statement cover the exact final buffers used by timing.
        VerifyParity(stripes, options, requested);
        const uint64_t checksum = ParityChecksum(stripes, options, requested);

        const uint64_t batch = static_cast<uint64_t>(options.batch);
        const uint64_t input_bytes = CheckedProduct(CheckedProduct(
            options.k, options.bytes, "logical input bytes"), batch,
            "logical input bytes");
        const uint64_t output_bytes = CheckedProduct(CheckedProduct(
            options.q, options.bytes, "requested parity bytes"), batch,
            "requested parity bytes");
        if (input_bytes > std::numeric_limits<uint64_t>::max() - output_bytes)
            Fail("logical input plus output bytes overflows uint64");
        const uint64_t direct_terms = CheckedProduct(CheckedProduct(
            options.k, options.q, "direct row terms"), batch,
            "direct row terms");
        const uint64_t direct_accumulations = CheckedProduct(CheckedProduct(
            options.k - 1, options.q, "direct accumulations"), batch,
            "direct accumulations");
        const uint64_t symbol_bytes =
            leo2_codec_field(codec) == LEO2_FIELD_GF16 ? 2 : 1;
        const uint64_t symbols_per_shard = options.bytes / symbol_bytes;
        const uint64_t direct_symbol_terms = CheckedProduct(
            direct_terms, symbols_per_shard, "direct field-symbol terms");
        const uint64_t direct_xor_symbols = CheckedProduct(
            direct_accumulations, symbols_per_shard,
            "direct accumulation symbols");
        const uint64_t direct_source_bytes = CheckedProduct(
            direct_terms, options.bytes, "modeled direct source bytes");
        const uint64_t direct_output_read_bytes = CheckedProduct(
            direct_accumulations, options.bytes,
            "modeled direct output read bytes");
        const uint64_t direct_output_write_bytes = CheckedProduct(
            direct_terms, options.bytes, "modeled direct output write bytes");
        const double amortized_us = encode.median_us +
            setup.median_us / static_cast<double>(options.reuse);

        std::ostringstream json;
        json.imbue(std::locale::classic());
        json << std::fixed << std::setprecision(6);
        json << "{\n"
             << "  \"schema\": \"leopard2-direct-encode-benchmark-v2\",\n"
             << "  \"build\": {\n"
             << "    \"compiler\": \"" << CompilerName() << "\",\n"
             << "    \"compiler_version\": \""
             << JsonEscape(CompilerVersion()) << "\",\n"
             << "    \"cplusplus\": " << __cplusplus << ",\n"
             << "    \"backend_variant\": \""
             << JsonEscape(LEO2_BENCHMARK_BUILD_VARIANT) << "\",\n"
             << "    \"build_type\": \""
             << JsonEscape(LEO2_BENCHMARK_BUILD_TYPE) << "\",\n"
             << "    \"build_configuration_sha256\": \""
             << LEO2_BENCHMARK_BUILD_CONFIGURATION_SHA256 << "\",\n"
             << "    \"test_hooks\": true,\n"
             << "    \"source_commit\": \""
             << LEO2_BENCHMARK_SOURCE_COMMIT << "\",\n"
             << "    \"source_tree\": \""
             << LEO2_BENCHMARK_SOURCE_TREE << "\",\n"
             << "    \"source_tracked_dirty\": "
             << (LEO2_BENCHMARK_SOURCE_TRACKED_DIRTY ? "true" : "false")
             << "\n"
             << "  },\n"
             << "  \"parameters\": {\n"
             << "    \"K\": " << options.k << ",\n"
             << "    \"R\": " << options.r << ",\n"
             << "    \"Q\": " << options.q << ",\n"
             << "    \"requested_parity_indices\": [";
        bool first = true;
        for (uint32_t i = 0; i < options.r; ++i)
        {
            if (!requested[i])
                continue;
            if (!first)
                json << ", ";
            json << i;
            first = false;
        }
        json << "],\n"
             << "    \"requested_profile\": \""
             << ProfileName(options.profile) << "\",\n"
             << "    \"requested_field\": \""
             << FieldName(options.field) << "\",\n"
             << "    \"forced_mode\": \"" << ModeName(options.mode) << "\",\n"
             << "    \"shard_bytes\": " << options.bytes << ",\n"
             << "    \"batch\": " << options.batch << ",\n"
             << "    \"reuse\": " << options.reuse << ",\n"
             << "    \"iterations\": " << options.iterations << ",\n"
             << "    \"warmups\": " << options.warmups << ",\n"
             << "    \"thread_count\": " << options.threads << ",\n"
             << "    \"seed\": " << options.seed << "\n"
             << "  },\n"
             << "  \"resolved\": {\n"
             << "    \"profile\": \""
             << ProfileName(leo2_codec_profile(codec)) << "\",\n"
             << "    \"field\": \""
             << FieldName(leo2_codec_field(codec)) << "\",\n"
             << "    \"backend\": \""
             << BackendName(leo2_context_backend(context)) << "\",\n"
             << "    \"thread_count\": "
             << leo2_context_thread_count(context) << ",\n"
             << "    \"parent_count\": "
             << leo2_codec_parent_count(codec) << ",\n"
             << "    \"padded_side\": "
             << leo2_codec_padded_side(codec) << ",\n"
             << "    \"direct_capable\": "
             << (direct_capable ? "true" : "false") << ",\n"
             << "    \"timed_path_is_direct\": "
             << (selected_direct ? "true" : "false") << "\n"
             << "  },\n"
             << "  \"correctness\": {\n"
             << "    \"selected_transform_reference_parity_match\": true,\n"
             << "    \"direct_transform_parity_match\": "
             << (selected_direct ? "true" : "null") << ",\n"
             << "    \"unrequested_outputs_untouched\": true,\n"
             << "    \"parity_checksum_fnv1a64\": \"0x"
             << std::hex << std::setw(16) << std::setfill('0') << checksum
             << std::dec << std::setfill(' ') << "\"\n"
             << "  },\n"
             << "  \"memory\": {\n"
             << "    \"scratch_alignment\": "
             << leo2_scratch_alignment() << ",\n"
             << "    \"encode_scratch_bytes_per_stripe\": "
             << scratch_bytes << ",\n"
             << "    \"reference_scratch_bytes_per_stripe\": "
             << reference_scratch_bytes << ",\n"
             << "    \"benchmark_allocated_scratch_bytes_per_stripe\": "
             << allocated_scratch_bytes << ",\n"
             << "    \"encode_scratch_bytes_batch\": "
             << CheckedProduct(scratch_bytes, batch, "batch scratch bytes")
             << ",\n"
             << "    \"benchmark_allocated_scratch_bytes_batch\": "
             << CheckedProduct(allocated_scratch_bytes, batch,
                    "benchmark allocated batch scratch bytes")
             << "\n"
             << "  },\n"
             << "  \"operation_model\": {\n"
             << "    \"direct_row_terms\": " << direct_terms << ",\n"
             << "    \"direct_output_initializations\": "
             << CheckedProduct(options.q, batch,
                    "direct output initializations") << ",\n"
             << "    \"direct_output_accumulations\": "
             << direct_accumulations << ",\n"
             << "    \"fixed_coefficient_symbol_terms_before_unit_"
             << "specialization\": " << direct_symbol_terms << ",\n"
             << "    \"xor_accumulation_symbols\": "
             << direct_xor_symbols << ",\n"
             << "    \"modeled_source_bytes_read\": "
             << direct_source_bytes << ",\n"
             << "    \"modeled_output_bytes_read\": "
             << direct_output_read_bytes << ",\n"
             << "    \"modeled_output_bytes_written\": "
             << direct_output_write_bytes << ",\n"
             << "    \"model_scope\": \"direct streaming kernels before "
             << "cache effects; unit coefficients specialize to copy/XOR\",\n"
             << "    \"direct_model_applies_to_timed_path\": "
             << (selected_direct ? "true" : "false") << ",\n"
             << "    \"transform_operation_counts\": null,\n"
             << "    \"hardware_counters\": null\n"
             << "  },\n"
             << "  \"metrics\": {\n"
             << "    \"codec_setup\": ";
        WriteSummary(json, setup, "", 4);
        json << ",\n    \"encode_execution\": ";
        WriteEncodeSummary(json, encode, input_bytes, output_bytes, 4);
        json << ",\n"
             << "    \"encode_amortized_at_reuse\": {\n"
             << "      \"reuse_count\": " << options.reuse << ",\n"
             << "      \"derived_median_us_per_batch_call\": "
             << amortized_us << ",\n"
             << "      \"logical_input_GB_per_s\": "
             << GigabytesPerSecond(input_bytes, amortized_us) << ",\n"
             << "      \"requested_parity_output_GB_per_s\": "
             << GigabytesPerSecond(output_bytes, amortized_us) << "\n"
             << "    }\n"
             << "  },\n"
             << "  \"methodology\": {\n"
             << "    \"codec_setup_scope\": "
             << "\"codec_create + encode-mode selection + codec_destroy; "
             << "context reused\",\n"
             << "    \"timed_encode_allocations\": "
             << "\"all benchmark-owned storage and sample vectors are "
             << "preallocated\",\n"
             << "    \"affinity\": \"pinning is handled externally\",\n"
             << "    \"rate_semantics\": "
             << "\"logical input counts K shards once; output counts only "
             << "the Q non-null requested parity shards\",\n"
             << "    \"counter_scope\": "
             << "\"direct analytic counts are reported above; transform "
             << "counts and optional perf counters remain in leopard-79h.16\",\n"
             << "    \"production_autotuning\": false\n"
             << "  }\n"
             << "}\n";

        if (options.output == "-")
            std::cout << json.str();
        else
        {
            std::ofstream output(
                options.output.c_str(), std::ios::out | std::ios::trunc);
            if (!output)
                Fail("cannot open JSON output: " + options.output);
            output << json.str();
            if (!output)
                Fail("failed writing JSON output: " + options.output);
        }

        leo2_codec_destroy(codec);
        codec = NULL;
        leo2_context_destroy(context);
        return 0;
    }
    catch (...)
    {
        leo2_codec_destroy(codec);
        leo2_context_destroy(context);
        throw;
    }
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
        std::cerr << "leopard2 direct encode benchmark: "
                  << error.what() << std::endl;
        return 1;
    }
}
