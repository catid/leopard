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

#include "leopard.h"
#include "leopard2.h"
#include "Leopard2Dispatch.h"
#if defined(LEO2_HIGH_LOW_DUALITY_ATTRIBUTION) || \
    defined(LEO2_HIGH_DECODE_COPY_ATTRIBUTION)
#include "Leopard2Direct.h"
#endif
#if defined(LEO2_BENCHMARK_SOURCE_ATTESTATION)
#if !defined(LEO2_BENCHMARK_SOURCE_ATTESTATION_HEADER)
#error "source-attested benchmark requires its exact generated header path"
#endif
#include LEO2_BENCHMARK_SOURCE_ATTESTATION_HEADER
#if !defined(LEO2_BENCHMARK_SOURCE_COMMIT) || \
    !defined(LEO2_BENCHMARK_SOURCE_TREE) || \
    !defined(LEO2_BENCHMARK_SOURCE_TRACKED_DIRTY)
#error "source-attested benchmark requires commit, tree, and dirty definitions"
#endif
#endif
#if defined(LEO2_HIGH_DECODE_COPY_ATTRIBUTION)
#include "LeopardFF8.h"
#include "LeopardFF16.h"
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

struct Options
{
    uint32_t k;
    uint32_t r;
    leo2_profile profile;
    leo2_field field;
    leo2_backend backend;
    uint64_t bytes;
    uint32_t losses;
    size_t batch;
    size_t reuse;
    size_t iterations;
    size_t warmup;
    uint32_t threads;
    uint64_t seed;
    bool force_generic_decode;
    bool force_specialized_decode;
    bool force_tiled_decode;
    bool force_materialized_decode;
    bool skip_legacy;
    bool retain_samples;
    bool report_decode_path;
    bool report_direct_executor;
    bool attest_source;
#if defined(LEO2_HIGH_LOW_DUALITY_ATTRIBUTION)
    bool force_translated_low;
    bool force_native_high;
#endif
#if defined(LEO2_HIGH_DECODE_COPY_ATTRIBUTION)
    enum HighEvaluatorMode
    {
        HighEvaluatorUnset,
        HighEvaluatorNoCopy,
        HighEvaluatorCopyFallback
    } high_evaluator_mode;
#endif
    std::string output;

    Options()
        : k(240)
        , r(16)
        , profile(LEO2_PROFILE_AUTO)
        , field(LEO2_FIELD_AUTO)
        , backend(LEO2_BACKEND_AUTO)
        , bytes(65536)
        , losses(1)
        , batch(1)
        , reuse(8)
        , iterations(9)
        , warmup(2)
        , threads(1)
        , seed(1)
        , force_generic_decode(false)
        , force_specialized_decode(false)
        , force_tiled_decode(false)
        , force_materialized_decode(false)
        , skip_legacy(false)
        , retain_samples(false)
        , report_decode_path(false)
        , report_direct_executor(false)
        , attest_source(false)
#if defined(LEO2_HIGH_LOW_DUALITY_ATTRIBUTION)
        , force_translated_low(false)
        , force_native_high(false)
#endif
#if defined(LEO2_HIGH_DECODE_COPY_ATTRIBUTION)
        , high_evaluator_mode(HighEvaluatorUnset)
#endif
        , output("-")
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
    AlignedBuffer()
        : data_(NULL)
        , size_(0)
    {}

    explicit AlignedBuffer(size_t bytes)
        : data_(NULL)
        , size_(0)
    {
        Reset(bytes);
    }

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

struct Stripe
{
    AlignedBuffer original_storage;
    AlignedBuffer recovery_storage;
    AlignedBuffer restored_storage;
    AlignedBuffer encode_scratch;
    AlignedBuffer decode_scratch;
    std::vector<const void*> original;
    std::vector<const void*> received_original;
    std::vector<void*> recovery_output;
    std::vector<const void*> received_recovery;
    std::vector<void*> restored;
};

struct LegacyStripe
{
    AlignedBuffer encode_storage;
    AlignedBuffer decode_storage;
    std::vector<void*> encode_work;
    std::vector<void*> decode_work;
    std::vector<const void*> received_original;
    std::vector<const void*> received_recovery;
};

struct Summary
{
    double median_us;
    double mad_us;
    double minimum_us;
    double maximum_us;
    std::vector<double> samples_us;
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

static size_t CheckedSize(uint64_t count, uint64_t bytes, const char* what)
{
    if (count != 0 && bytes > std::numeric_limits<size_t>::max() / count)
        Fail(std::string(what) + " size overflow");
    return static_cast<size_t>(count * bytes);
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
    while (digit_count < text.size() && text[digit_count] >= '0' && text[digit_count] <= '9')
        ++digit_count;
    if (digit_count == 0)
        Fail("invalid value for --bytes: " + text);
    const uint64_t base = ParseUnsigned(text.substr(0, digit_count), "--bytes");
    std::string suffix = text.substr(digit_count);
    std::transform(suffix.begin(), suffix.end(), suffix.begin(),
        [](unsigned char character) { return static_cast<char>(std::tolower(character)); });
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
        multiplier = UINT64_C(1024) * 1024;
    else if (suffix == "gi" || suffix == "gib")
        multiplier = UINT64_C(1024) * 1024 * 1024;
    else
        Fail("invalid suffix for --bytes: " + suffix);
    if (base > std::numeric_limits<uint64_t>::max() / multiplier)
        Fail("--bytes overflows uint64");
    return base * multiplier;
}

static leo2_profile ParseProfile(const std::string& text)
{
    if (text == "auto") return LEO2_PROFILE_AUTO;
    if (text == "high" || text == "legacy-high") return LEO2_PROFILE_LEGACY_HIGH_V1;
    if (text == "low") return LEO2_PROFILE_LOW_V1;
    if (text == "exact") return LEO2_PROFILE_EXACT_EXPERIMENTAL_V1;
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

static leo2_backend ParseBackend(const std::string& text)
{
    if (text == "auto") return LEO2_BACKEND_AUTO;
    if (text == "scalar") return LEO2_BACKEND_SCALAR;
    if (text == "ssse3") return LEO2_BACKEND_SSSE3;
    if (text == "avx2") return LEO2_BACKEND_AVX2;
    if (text == "avx512" || text == "avx512vl")
        return LEO2_BACKEND_AVX512;
    if (text == "gfni" || text == "avx2-gfni")
        return LEO2_BACKEND_GFNI;
    if (text == "neon") return LEO2_BACKEND_NEON;
    Fail("invalid --backend: " + text);
    return LEO2_BACKEND_AUTO;
}

static const char* NeedValue(int argc, char** argv, int& index)
{
    if (++index >= argc)
        Fail(std::string("missing value after ") + argv[index - 1]);
    return argv[index];
}

static void Usage(std::ostream& output, const char* program)
{
    output
        << "Usage: " << program << " [options]\n"
        << "  --k N                 Original shard count (default 240)\n"
        << "  --r N                 Recovery shard count (default 16)\n"
        << "  --profile NAME        auto, high, low, exact (default auto)\n"
        << "  --field NAME          auto, gf8, gf16 (default auto)\n"
        << "  --backend NAME        auto, scalar, ssse3, avx2, avx512, gfni,\n"
        << "                        neon\n"
        << "  --bytes N[KiB|MiB]    Bytes per shard (default 64KiB)\n"
        << "  --loss N              Missing original shards (default 1)\n"
        << "  --batch N             Independent stripes per call (default 1)\n"
        << "  --reuse N             Calls per timing sample (default 8)\n"
        << "  --iterations N        Timing samples (default 9)\n"
        << "  --warmup N            Untimed calls (default 2)\n"
        << "  --threads N           Context thread count (default 1)\n"
        << "  --seed N              Deterministic seed (default 1)\n"
        << "  --force-generic       Use the retained O(N log N) decoder\n"
        << "  --force-specialized   Use the profile-specific transform decoder\n"
        << "  --force-tiled         Use the side-sized specialized kernel\n"
        << "  --force-materialized  Use the retained N-slot specialized kernel\n"
        << "  --skip-legacy         Do not run the in-tree legacy comparison\n"
        << "  --retain-samples      Emit raw timing samples using benchmark schema v2\n"
        << "  --report-decode-path  Emit internal selected-path metadata using schema v3\n"
#if !defined(LEO2_HIGH_DECODE_COPY_ATTRIBUTION) && \
    !defined(LEO2_HIGH_LOW_DUALITY_ATTRIBUTION)
        << "  --report-direct-executor\n"
        << "                         Add direct loop-order metadata using schema v6\n"
#endif
#if defined(LEO2_BENCHMARK_SOURCE_ATTESTATION)
        << "  --attest-source       Embed committed source identity using schema v5\n"
#endif
#if defined(LEO2_HIGH_DECODE_COPY_ATTRIBUTION)
        << "  --high-evaluator-mode NAME\n"
        << "                         Attribution-only: no-copy or copy-fallback\n"
#endif
#if defined(LEO2_HIGH_LOW_DUALITY_ATTRIBUTION)
        << "  --force-translated-low\n"
        << "                         Attribution-only: run high wire through Algorithm 4\n"
        << "  --force-native-high\n"
        << "                         Attribution-only: run high wire through Algorithm 5\n"
#endif
        << "  --json PATH           JSON output path, or - for stdout\n"
        << "  --help                 Show this message\n";
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
        else if (argument == "--k") options.k = ParseUint32(NeedValue(argc, argv, i), "--k");
        else if (argument == "--r") options.r = ParseUint32(NeedValue(argc, argv, i), "--r");
        else if (argument == "--profile") options.profile = ParseProfile(NeedValue(argc, argv, i));
        else if (argument == "--field") options.field = ParseField(NeedValue(argc, argv, i));
        else if (argument == "--backend") options.backend = ParseBackend(NeedValue(argc, argv, i));
        else if (argument == "--bytes") options.bytes = ParseBytes(NeedValue(argc, argv, i));
        else if (argument == "--loss" || argument == "--losses") options.losses = ParseUint32(NeedValue(argc, argv, i), "--loss");
        else if (argument == "--batch") options.batch = ParseSize(NeedValue(argc, argv, i), "--batch");
        else if (argument == "--reuse") options.reuse = ParseSize(NeedValue(argc, argv, i), "--reuse");
        else if (argument == "--iterations") options.iterations = ParseSize(NeedValue(argc, argv, i), "--iterations");
        else if (argument == "--warmup") options.warmup = ParseSize(NeedValue(argc, argv, i), "--warmup");
        else if (argument == "--threads" || argument == "--thread-count") options.threads = ParseUint32(NeedValue(argc, argv, i), "--threads");
        else if (argument == "--seed") options.seed = ParseUnsigned(NeedValue(argc, argv, i), "--seed");
        else if (argument == "--force-generic") options.force_generic_decode = true;
        else if (argument == "--force-specialized") options.force_specialized_decode = true;
        else if (argument == "--force-tiled") options.force_tiled_decode = true;
        else if (argument == "--force-materialized") options.force_materialized_decode = true;
        else if (argument == "--skip-legacy") options.skip_legacy = true;
        else if (argument == "--retain-samples") options.retain_samples = true;
        else if (argument == "--report-decode-path") options.report_decode_path = true;
        else if (argument == "--report-direct-executor")
        {
#if defined(LEO2_HIGH_DECODE_COPY_ATTRIBUTION) || \
    defined(LEO2_HIGH_LOW_DUALITY_ATTRIBUTION)
            Fail("--report-direct-executor is unavailable in the "
                 "attribution benchmark");
#else
            options.report_direct_executor = true;
            options.report_decode_path = true;
#endif
        }
#if defined(LEO2_BENCHMARK_SOURCE_ATTESTATION)
        else if (argument == "--attest-source") options.attest_source = true;
#endif
#if defined(LEO2_HIGH_DECODE_COPY_ATTRIBUTION)
        else if (argument == "--high-evaluator-mode")
        {
            const std::string mode = NeedValue(argc, argv, i);
            if (mode == "no-copy")
                options.high_evaluator_mode = Options::HighEvaluatorNoCopy;
            else if (mode == "copy-fallback")
                options.high_evaluator_mode = Options::HighEvaluatorCopyFallback;
            else
                Fail("invalid --high-evaluator-mode: " + mode);
        }
#endif
#if defined(LEO2_HIGH_LOW_DUALITY_ATTRIBUTION)
        else if (argument == "--force-translated-low")
            options.force_translated_low = true;
        else if (argument == "--force-native-high")
            options.force_native_high = true;
#endif
        else if (argument == "--json" || argument == "--output") options.output = NeedValue(argc, argv, i);
        else Fail("unknown argument: " + argument);
    }

    if (options.k == 0 || options.r == 0)
        Fail("--k and --r must be positive and fit in uint32");
    if (options.bytes == 0 || options.bytes > std::numeric_limits<size_t>::max())
        Fail("--bytes must be positive and fit in size_t");
    if (options.batch == 0 || options.reuse == 0 || options.iterations == 0 || options.threads == 0)
        Fail("--batch, --reuse, --iterations, and --threads must be positive");
    if (options.losses > options.k)
        Fail("--loss cannot exceed K");
    if (options.losses > options.r)
        Fail("--loss cannot exceed R when only transmitted recovery shards are used");
    if (options.force_generic_decode && options.force_specialized_decode)
        Fail("--force-generic and --force-specialized are mutually exclusive");
    if (options.force_tiled_decode && options.force_materialized_decode)
        Fail("--force-tiled and --force-materialized are mutually exclusive");
    if (options.force_generic_decode &&
        (options.force_tiled_decode || options.force_materialized_decode))
        Fail("--force-generic cannot select a specialized workspace kernel");
#if defined(LEO2_BENCHMARK_SOURCE_ATTESTATION)
    if (options.attest_source && options.report_direct_executor)
        Fail("--attest-source and --report-direct-executor use distinct "
             "JSON schemas");
    if (options.attest_source && options.report_decode_path)
        Fail("--attest-source and --report-decode-path use distinct JSON schemas");
#endif
#if defined(LEO2_HIGH_DECODE_COPY_ATTRIBUTION)
    if (options.attest_source)
        Fail("the attribution benchmark does not support --attest-source");
    if (options.high_evaluator_mode == Options::HighEvaluatorUnset)
        Fail("the attribution benchmark requires --high-evaluator-mode");
    if (options.profile != LEO2_PROFILE_LEGACY_HIGH_V1 ||
        (options.field != LEO2_FIELD_GF8 && options.field != LEO2_FIELD_GF16) ||
        !options.force_specialized_decode || options.force_generic_decode ||
        options.force_tiled_decode == options.force_materialized_decode ||
        !options.skip_legacy || !options.retain_samples ||
        !options.report_decode_path || options.threads != 1)
        Fail("the attribution benchmark requires explicit high profile/field, "
             "one forced specialized workspace, --skip-legacy, "
             "--retain-samples, --report-decode-path, and one thread");
#endif
#if defined(LEO2_HIGH_LOW_DUALITY_ATTRIBUTION)
    if (options.attest_source)
        Fail("the high/low duality benchmark does not support --attest-source");
    if (options.force_translated_low && options.force_native_high)
        Fail("translated-low and native-high controls are mutually exclusive");
    if ((options.force_translated_low || options.force_native_high) &&
        (options.profile != LEO2_PROFILE_LEGACY_HIGH_V1 ||
         (options.field != LEO2_FIELD_GF8 &&
          options.field != LEO2_FIELD_GF16) ||
         options.force_generic_decode ||
         options.force_tiled_decode == options.force_materialized_decode ||
         !options.retain_samples || !options.report_decode_path ||
         options.threads != 1))
    {
        Fail("translated-low attribution requires explicit high profile/field, "
             "one forced workspace, --retain-samples, "
             "--report-decode-path, and one thread");
    }
#endif
    return options;
}

#if defined(LEO2_HIGH_DECODE_COPY_ATTRIBUTION)
static const char* HighEvaluatorModeName(Options::HighEvaluatorMode mode)
{
    switch (mode)
    {
    case Options::HighEvaluatorNoCopy: return "no-copy";
    case Options::HighEvaluatorCopyFallback: return "copy-fallback";
    case Options::HighEvaluatorUnset: break;
    }
    return "unset";
}

struct HighEvaluatorAttribution
{
    uint64_t output_blocks;
    uint64_t butterfly2_out_of_place;
    uint64_t butterfly4_out_of_place;
    uint64_t compatibility_copy_fallbacks;
    uint64_t pruned_output_blocks;
    uint64_t mature_output_blocks;
};

static void PrepareHighEvaluatorAttribution(const Options& options)
{
    const bool copy = options.high_evaluator_mode ==
        Options::HighEvaluatorCopyFallback;
    leopard::ff8::TestOnlySetHighDecodeCopyFallback(copy);
    leopard::ff16::TestOnlySetHighDecodeCopyFallback(copy);
    leopard::ff8::TestOnlyResetHighDecodeCounts();
    leopard::ff16::TestOnlyResetHighDecodeCounts();
    if (leopard::ff8::TestOnlyHighDecodeCopyFallbackEnabled() != copy ||
        leopard::ff16::TestOnlyHighDecodeCopyFallbackEnabled() != copy)
        Fail("high evaluator attribution selector did not latch");
}

static HighEvaluatorAttribution ReadHighEvaluatorAttribution(
    leo2_field field,
    const Options& options)
{
    const leopard::ff8::TestOnlyHighDecodeCounts ff8 =
        leopard::ff8::TestOnlyGetHighDecodeCounts();
    const leopard::ff16::TestOnlyHighDecodeCounts ff16 =
        leopard::ff16::TestOnlyGetHighDecodeCounts();
    const bool selected_ff8 = field == LEO2_FIELD_GF8;
    const uint64_t inactive_output_blocks = selected_ff8
        ? ff16.output_blocks : ff8.output_blocks;
    const uint64_t inactive_copy_fallbacks = selected_ff8
        ? ff16.compatibility_copy_fallbacks
        : ff8.compatibility_copy_fallbacks;
    const uint64_t inactive_out_of_place = selected_ff8
        ? ff16.fft_butterfly2_out_of_place +
            ff16.fft_butterfly4_out_of_place
        : ff8.fft_butterfly2_out_of_place +
            ff8.fft_butterfly4_out_of_place;
    const uint64_t inactive_classified_blocks = selected_ff8
        ? ff16.pruned_output_blocks + ff16.mature_output_blocks
        : ff8.pruned_output_blocks + ff8.mature_output_blocks;
    if (inactive_output_blocks != 0 || inactive_copy_fallbacks != 0 ||
        inactive_out_of_place != 0 || inactive_classified_blocks != 0)
        Fail("inactive field recorded high evaluator work");

    HighEvaluatorAttribution result;
    if (selected_ff8)
    {
        result.output_blocks = ff8.output_blocks;
        result.butterfly2_out_of_place = ff8.fft_butterfly2_out_of_place;
        result.butterfly4_out_of_place = ff8.fft_butterfly4_out_of_place;
        result.compatibility_copy_fallbacks =
            ff8.compatibility_copy_fallbacks;
        result.pruned_output_blocks = ff8.pruned_output_blocks;
        result.mature_output_blocks = ff8.mature_output_blocks;
    }
    else
    {
        result.output_blocks = ff16.output_blocks;
        result.butterfly2_out_of_place = ff16.fft_butterfly2_out_of_place;
        result.butterfly4_out_of_place = ff16.fft_butterfly4_out_of_place;
        result.compatibility_copy_fallbacks =
            ff16.compatibility_copy_fallbacks;
        result.pruned_output_blocks = ff16.pruned_output_blocks;
        result.mature_output_blocks = ff16.mature_output_blocks;
    }
    const uint64_t out_of_place = result.butterfly2_out_of_place +
        result.butterfly4_out_of_place;
    if (result.output_blocks == 0)
        Fail("attribution benchmark executed no Algorithm 5 output blocks");
    if (result.pruned_output_blocks + result.mature_output_blocks !=
        result.output_blocks)
        Fail("Algorithm 5 output-block classification is incomplete");
    if (options.high_evaluator_mode == Options::HighEvaluatorNoCopy)
    {
        if (result.compatibility_copy_fallbacks != 0 || out_of_place == 0)
            Fail("no-copy attribution did not exclusively use out-of-place evaluation");
    }
    else if (result.compatibility_copy_fallbacks != result.output_blocks ||
             out_of_place != 0)
        Fail("copy attribution did not force one fallback per output block");
    return result;
}
#endif

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

static std::string JsonEscape(const std::string& input)
{
    std::ostringstream output;
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
                output << "\\u" << std::hex << std::setw(4) << std::setfill('0')
                       << static_cast<unsigned>(character) << std::dec;
            else
                output << input[i];
        }
    }
    return output.str();
}

static double Median(std::vector<double> values)
{
    const size_t middle = values.size() / 2;
    std::nth_element(values.begin(), values.begin() + middle, values.end());
    const double upper = values[middle];
    if ((values.size() & 1u) != 0)
        return upper;
    std::nth_element(values.begin(), values.begin() + middle - 1, values.begin() + middle);
    return (values[middle - 1] + upper) * 0.5;
}

static Summary Summarize(
    const std::vector<double>& samples,
    bool retain_samples)
{
    Summary summary;
    if (retain_samples)
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
static Summary Measure(
    size_t iterations,
    size_t inner_calls,
    bool retain_samples,
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
    return Summarize(samples, retain_samples);
}

static void FillOriginals(Stripe& stripe, const Options& options, size_t stripe_index)
{
    XorShift64 random(options.seed ^
        (UINT64_C(0x9e3779b97f4a7c15) * (stripe_index + 1)));
    uint8_t* output = stripe.original_storage.bytes();
    const size_t total = CheckedSize(options.k, options.bytes, "original data");
    for (size_t i = 0; i < total; ++i)
        output[i] = static_cast<uint8_t>(random.Next() >> 56);
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

static void InitializeStripe(
    Stripe& stripe,
    const Options& options,
    const std::vector<uint32_t>& losses,
    size_t stripe_index,
    size_t encode_scratch_bytes,
    size_t decode_scratch_bytes)
{
    stripe.original_storage.Reset(CheckedSize(options.k, options.bytes, "original data"));
    stripe.recovery_storage.Reset(CheckedSize(options.r, options.bytes, "recovery data"));
    stripe.restored_storage.Reset(CheckedSize(options.k, options.bytes, "restored data"));
    stripe.encode_scratch.Reset(encode_scratch_bytes);
    stripe.decode_scratch.Reset(decode_scratch_bytes);
    stripe.original.resize(options.k);
    stripe.received_original.resize(options.k);
    stripe.recovery_output.resize(options.r);
    stripe.received_recovery.resize(options.r);
    stripe.restored.assign(options.k, static_cast<void*>(NULL));
    FillOriginals(stripe, options, stripe_index);

    for (uint32_t i = 0; i < options.k; ++i)
    {
        uint8_t* data = stripe.original_storage.bytes() +
            static_cast<size_t>(i) * static_cast<size_t>(options.bytes);
        stripe.original[i] = data;
        stripe.received_original[i] = IsLost(losses, i) ? NULL : data;
        if (IsLost(losses, i))
            stripe.restored[i] = stripe.restored_storage.bytes() +
                static_cast<size_t>(i) * static_cast<size_t>(options.bytes);
    }
    for (uint32_t i = 0; i < options.r; ++i)
    {
        uint8_t* data = stripe.recovery_storage.bytes() +
            static_cast<size_t>(i) * static_cast<size_t>(options.bytes);
        stripe.recovery_output[i] = data;
        stripe.received_recovery[i] = data;
    }
}

static void CheckRestored(
    const Stripe& stripe,
    const Options& options,
    const std::vector<uint32_t>& losses,
    const char* implementation)
{
    for (size_t loss_i = 0; loss_i < losses.size(); ++loss_i)
    {
        const uint32_t index = losses[loss_i];
        const uint8_t* expected = stripe.original_storage.bytes() +
            static_cast<size_t>(index) * static_cast<size_t>(options.bytes);
        const uint8_t* actual = stripe.restored_storage.bytes() +
            static_cast<size_t>(index) * static_cast<size_t>(options.bytes);
        if (memcmp(expected, actual, static_cast<size_t>(options.bytes)) != 0)
        {
            std::ostringstream message;
            message << implementation << " restored shard " << index << " incorrectly";
            Fail(message.str());
        }
    }
}

static void InitializeLegacyStripe(
    LegacyStripe& legacy,
    const Stripe& stripe,
    const Options& options,
    unsigned encode_work_count,
    unsigned decode_work_count)
{
    legacy.encode_storage.Reset(CheckedSize(encode_work_count, options.bytes, "legacy encode work"));
    legacy.decode_storage.Reset(CheckedSize(decode_work_count, options.bytes, "legacy decode work"));
    legacy.encode_work.resize(encode_work_count);
    legacy.decode_work.resize(decode_work_count);
    legacy.received_original = stripe.received_original;
    legacy.received_recovery.resize(options.r);
    for (unsigned i = 0; i < encode_work_count; ++i)
        legacy.encode_work[i] = legacy.encode_storage.bytes() +
            static_cast<size_t>(i) * static_cast<size_t>(options.bytes);
    for (unsigned i = 0; i < decode_work_count; ++i)
        legacy.decode_work[i] = legacy.decode_storage.bytes() +
            static_cast<size_t>(i) * static_cast<size_t>(options.bytes);
    for (uint32_t i = 0; i < options.r; ++i)
        legacy.received_recovery[i] = legacy.encode_work[i];
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

static double GigabytesPerSecond(uint64_t bytes, double microseconds)
{
    if (bytes == 0 || microseconds <= 0)
        return 0;
    return static_cast<double>(bytes) / (microseconds * 1000.0);
}

static uint64_t Fnv1a64Update(uint64_t digest, const void* data, size_t bytes)
{
    static const uint64_t kPrime = UINT64_C(1099511628211);
    const uint8_t* input = static_cast<const uint8_t*>(data);
    for (size_t i = 0; i < bytes; ++i)
    {
        digest ^= input[i];
        digest *= kPrime;
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

static uint64_t CheckedU64Product(uint64_t left, uint64_t right, const char* what)
{
    if (left != 0 && right > std::numeric_limits<uint64_t>::max() / left)
        Fail(std::string(what) + " overflows uint64");
    return left * right;
}

static uint64_t CheckedU64Sum(uint64_t left, uint64_t right, const char* what)
{
    if (right > std::numeric_limits<uint64_t>::max() - left)
        Fail(std::string(what) + " overflows uint64");
    return left + right;
}

static void WriteOptionalRate(std::ostream& output, uint64_t bytes, double microseconds)
{
    if (bytes == 0)
        output << "null";
    else
        output << GigabytesPerSecond(bytes, microseconds);
}

static void WriteSummary(
    std::ostream& output,
    const Summary& summary,
    uint64_t input_bytes,
    const char* input_name,
    uint64_t output_bytes,
    const char* output_name,
    unsigned indent,
    bool retain_samples)
{
    const std::string spaces(indent, ' ');
    output << "{\n"
           << spaces << "  \"median_us_per_batch_call\": " << summary.median_us << ",\n"
           << spaces << "  \"mad_us_per_batch_call\": " << summary.mad_us << ",\n"
           << spaces << "  \"minimum_us_per_batch_call\": " << summary.minimum_us << ",\n"
           << spaces << "  \"maximum_us_per_batch_call\": " << summary.maximum_us << ",\n";
    if (retain_samples)
    {
        output << spaces << "  \"samples_us_per_batch_call\": [";
        for (size_t i = 0; i < summary.samples_us.size(); ++i)
        {
            if (i != 0)
                output << ", ";
            output << summary.samples_us[i];
        }
        output << "],\n";
    }
    output << spaces << "  \"" << input_name << "\": ";
    WriteOptionalRate(output, input_bytes, summary.median_us);
    output << ",\n" << spaces << "  \"" << output_name << "\": ";
    WriteOptionalRate(output, output_bytes, summary.median_us);
    output << "\n" << spaces << '}';
}

static void WriteSetupSummary(
    std::ostream& output,
    const Summary& summary,
    unsigned indent,
    bool retain_samples)
{
    const std::string spaces(indent, ' ');
    output << "{\n"
           << spaces << "  \"median_us\": " << summary.median_us << ",\n"
           << spaces << "  \"mad_us\": " << summary.mad_us << ",\n"
           << spaces << "  \"minimum_us\": " << summary.minimum_us << ",\n"
           << spaces << "  \"maximum_us\": " << summary.maximum_us;
    if (retain_samples)
    {
        output << ",\n" << spaces << "  \"samples_us\": [";
        for (size_t i = 0; i < summary.samples_us.size(); ++i)
        {
            if (i != 0)
                output << ", ";
            output << summary.samples_us[i];
        }
        output << ']';
    }
    output << "\n" << spaces << '}';
}

static void WriteAmortizedDecodeSummary(
    std::ostream& output,
    const Summary& setup,
    const Summary& execution,
    size_t reuse,
    uint64_t offered_input_bytes,
    uint64_t repaired_output_bytes,
    unsigned indent)
{
    const std::string spaces(indent, ' ');
    const double microseconds = execution.median_us +
        setup.median_us / static_cast<double>(reuse);
    output << "{\n"
           << spaces << "  \"reuse_count\": " << reuse << ",\n"
           << spaces << "  \"derived_median_us_per_batch_call\": " << microseconds << ",\n"
           << spaces << "  \"offered_received_GB_per_s\": ";
    WriteOptionalRate(output, offered_input_bytes, microseconds);
    output << ",\n" << spaces << "  \"repaired_output_GB_per_s\": ";
    WriteOptionalRate(output, repaired_output_bytes, microseconds);
    output << "\n" << spaces << '}';
}

static std::string LegacyUnavailableReason(
    const Options& options,
    const leo2_codec* codec)
{
    if (leo2_codec_profile(codec) != LEO2_PROFILE_LEGACY_HIGH_V1)
        return "old Leopard only defines the legacy high wire profile";
    if ((options.bytes & 63u) != 0)
        return "old Leopard requires shard bytes divisible by 64";
    if (options.r > options.k)
        return "old Leopard requires R <= K";
    const leo2_field legacy_field = leo2_codec_parent_count(codec) <= 256
        ? LEO2_FIELD_GF8 : LEO2_FIELD_GF16;
    if (leo2_codec_field(codec) != legacy_field)
        return "the selected field differs from old Leopard's automatic field";
    if (leo_encode_work_count(options.k, options.r) == 0 ||
        leo_decode_work_count(options.k, options.r) == 0)
        return "old Leopard rejected the requested counts";
    return std::string();
}

static int Run(const Options& options)
{
    leo2_context_options context_options;
    memset(&context_options, 0, sizeof(context_options));
    context_options.struct_size = sizeof(context_options);
    context_options.backend = options.backend;
    context_options.thread_count = options.threads;
    leo2_context* context = NULL;
    RequireLeo2(leo2_context_create(&context_options, &context), "context create");

    leo2_codec* codec = NULL;
    leo2_codec_options codec_options;
    memset(&codec_options, 0, sizeof(codec_options));
    codec_options.struct_size = sizeof(codec_options);
    codec_options.flags = 0;
    if (options.force_generic_decode)
        codec_options.flags |= LEO2_CODEC_FORCE_GENERIC_DECODE;
    if (options.force_specialized_decode)
        codec_options.flags |= LEO2_CODEC_FORCE_SPECIALIZED_DECODE;
    if (options.force_tiled_decode)
        codec_options.flags |= LEO2_CODEC_FORCE_TILED_DECODE;
    if (options.force_materialized_decode)
        codec_options.flags |= LEO2_CODEC_FORCE_MATERIALIZED_DECODE;
    RequireLeo2(leo2_codec_create(
        context, options.k, options.r, options.profile, options.field,
        &codec_options, &codec),
        "codec create");
#if defined(LEO2_HIGH_LOW_DUALITY_ATTRIBUTION)
    if (options.force_translated_low || options.force_native_high)
    {
        RequireLeo2(leo2_test_codec_set_decode_mode(
            codec, options.force_translated_low
                ? LEO2_TEST_DECODE_FORCE_TRANSLATED_LOW
                : LEO2_TEST_DECODE_FORCE_NATIVE_HIGH),
            "force attributed decoder");
    }
#endif
#if defined(LEO2_HIGH_DECODE_COPY_ATTRIBUTION)
    if (leo2_test_codec_translated_low_capable(codec))
    {
        RequireLeo2(leo2_test_codec_set_decode_mode(
            codec, LEO2_TEST_DECODE_FORCE_NATIVE_HIGH),
            "force native high evaluator attribution");
    }
#endif

    size_t encode_scratch_bytes = 0;
    size_t decode_scratch_bytes = 0;
    RequireLeo2(leo2_encode_scratch_size(codec, options.bytes, &encode_scratch_bytes),
        "encode scratch query");

    const std::vector<uint32_t> losses = SelectLosses(options);
    std::vector<uint8_t> original_present(options.k, 1);
    std::vector<uint8_t> recovery_present(options.r, 1);
    for (size_t i = 0; i < losses.size(); ++i)
        original_present[losses[i]] = 0;

    leo2_decode_plan* plan = NULL;
    RequireLeo2(leo2_decode_plan_create(
        codec, &original_present[0], &recovery_present[0], &plan), "decode plan create");
#if defined(LEO2_HIGH_LOW_DUALITY_ATTRIBUTION)
    if ((options.force_translated_low || options.force_native_high) &&
        leo2_test_decode_plan_uses_translated_low(plan) !=
            (options.force_translated_low ? 1 : 0))
        Fail("decode plan did not capture translated-low attribution mode");
#endif
    RequireLeo2(leo2_decode_plan_scratch_size(plan, options.bytes, &decode_scratch_bytes),
        "decode scratch query");

    leopard2_internal::DecodePathInfo decode_path_info;
    if (options.report_decode_path)
    {
        RequireLeo2(leopard2_internal::GetDecodePlanPathInfo(
            plan, options.bytes, options.batch > 1, &decode_path_info),
            "decode path introspection");
    }

#if defined(LEO2_HIGH_DECODE_COPY_ATTRIBUTION)
    PrepareHighEvaluatorAttribution(options);
#endif

    std::vector<std::unique_ptr<Stripe> > stripes;
    stripes.reserve(options.batch);
    std::vector<leo2_encode_batch_item> encode_items(options.batch);
    std::vector<leo2_decode_batch_item> decode_items(options.batch);
    for (size_t i = 0; i < options.batch; ++i)
    {
        stripes.push_back(std::unique_ptr<Stripe>(new Stripe));
        InitializeStripe(*stripes.back(), options, losses, i,
            encode_scratch_bytes, decode_scratch_bytes);
        Stripe& stripe = *stripes.back();
        encode_items[i].shard_bytes = options.bytes;
        encode_items[i].original = &stripe.original[0];
        encode_items[i].recovery = &stripe.recovery_output[0];
        encode_items[i].scratch = stripe.encode_scratch.data();
        encode_items[i].scratch_bytes = stripe.encode_scratch.size();
        decode_items[i].shard_bytes = options.bytes;
        decode_items[i].original = &stripe.received_original[0];
        decode_items[i].recovery = &stripe.received_recovery[0];
        decode_items[i].restored_original = &stripe.restored[0];
        decode_items[i].scratch = stripe.decode_scratch.data();
        decode_items[i].scratch_bytes = stripe.decode_scratch.size();
    }

    size_t encode_batch_preflight_bytes = 0;
    size_t decode_batch_preflight_bytes = 0;
    RequireLeo2(leo2_encode_batch_preflight_scratch_size(
        codec, encode_items.size(), &encode_batch_preflight_bytes),
        "encode batch preflight scratch query");
    RequireLeo2(leo2_decode_plan_batch_preflight_scratch_size(
        plan, decode_items.size(), &decode_batch_preflight_bytes),
        "decode batch preflight scratch query");
    AlignedBuffer encode_batch_preflight(encode_batch_preflight_bytes);
    AlignedBuffer decode_batch_preflight(decode_batch_preflight_bytes);
    const auto run_encode_batch = [&]() {
        const leo2_result result = encode_batch_preflight_bytes == 0
            ? leo2_encode_batch(
                codec, &encode_items[0], encode_items.size())
            : leo2_encode_batch_with_preflight_scratch(
                codec, &encode_items[0], encode_items.size(),
                encode_batch_preflight.data(),
                encode_batch_preflight.size());
        RequireLeo2(result, "encode batch");
    };
    const auto run_decode_batch = [&]() {
        const leo2_result result = decode_batch_preflight_bytes == 0
            ? leo2_decode_plan_execute_batch(
                plan, &decode_items[0], decode_items.size())
            : leo2_decode_plan_execute_batch_with_preflight_scratch(
                plan, &decode_items[0], decode_items.size(),
                decode_batch_preflight.data(),
                decode_batch_preflight.size());
        RequireLeo2(result, "decode batch");
    };

    const bool extended_schema = options.skip_legacy || options.retain_samples ||
        options.report_decode_path || options.report_direct_executor ||
        options.attest_source;
    const unsigned schema_version =
#if defined(LEO2_HIGH_DECODE_COPY_ATTRIBUTION)
        4;
#else
        options.report_direct_executor ? 6 :
        options.attest_source ? 5 :
        (options.report_decode_path ? 3 : (extended_schema ? 2 : 1));
#endif
    static const uint64_t kFnv1a64Offset = UINT64_C(14695981039346656037);
    uint64_t original_digest = kFnv1a64Offset;
    uint64_t parity_digest = kFnv1a64Offset;
    uint64_t recovered_digest = kFnv1a64Offset;

    run_encode_batch();
    if (extended_schema)
    {
        for (size_t stripe_index = 0; stripe_index < stripes.size(); ++stripe_index)
        {
            const Stripe& stripe = *stripes[stripe_index];
            original_digest = Fnv1a64Update(original_digest,
                stripe.original_storage.data(),
                CheckedSize(options.k, options.bytes, "original digest"));
            parity_digest = Fnv1a64Update(parity_digest,
                stripe.recovery_storage.data(),
                CheckedSize(options.r, options.bytes, "parity digest"));
        }
    }

    run_decode_batch();
    for (size_t i = 0; i < stripes.size(); ++i)
        CheckRestored(*stripes[i], options, losses, "Leopard2");
    if (extended_schema)
    {
        for (size_t stripe_index = 0; stripe_index < stripes.size(); ++stripe_index)
        {
            const Stripe& stripe = *stripes[stripe_index];
            for (size_t loss_i = 0; loss_i < losses.size(); ++loss_i)
            {
                const uint32_t index = losses[loss_i];
                const uint8_t* restored = stripe.restored_storage.bytes() +
                    static_cast<size_t>(index) * static_cast<size_t>(options.bytes);
                recovered_digest = Fnv1a64Update(recovered_digest, restored,
                    static_cast<size_t>(options.bytes));
            }
        }
    }

    const Summary codec_setup = Measure(
        options.iterations, 1, options.retain_samples, [&]() {
        leo2_codec* temporary = NULL;
        RequireLeo2(leo2_codec_create(
            context, options.k, options.r, options.profile, options.field,
            &codec_options, &temporary),
            "timed codec create");
#if defined(LEO2_HIGH_LOW_DUALITY_ATTRIBUTION)
        if (options.force_translated_low || options.force_native_high)
        {
            RequireLeo2(leo2_test_codec_set_decode_mode(
                temporary, options.force_translated_low
                    ? LEO2_TEST_DECODE_FORCE_TRANSLATED_LOW
                    : LEO2_TEST_DECODE_FORCE_NATIVE_HIGH),
                "timed force attributed decoder");
        }
#endif
#if defined(LEO2_HIGH_DECODE_COPY_ATTRIBUTION)
        if (leo2_test_codec_translated_low_capable(temporary))
        {
            RequireLeo2(leo2_test_codec_set_decode_mode(
                temporary, LEO2_TEST_DECODE_FORCE_NATIVE_HIGH),
                "timed force native high evaluator attribution");
        }
#endif
        leo2_codec_destroy(temporary);
    });

    const Summary plan_setup = Measure(
        options.iterations, 1, options.retain_samples, [&]() {
        leo2_decode_plan* temporary = NULL;
        RequireLeo2(leo2_decode_plan_create(codec,
            &original_present[0], &recovery_present[0], &temporary),
            "timed decode plan create");
        leo2_decode_plan_destroy(temporary);
    });

    for (size_t i = 0; i < options.warmup; ++i)
    {
        run_encode_batch();
        run_decode_batch();
    }

    const Summary encode_execution = Measure(
        options.iterations, options.reuse, options.retain_samples, [&]() {
        run_encode_batch();
    });
    const Summary decode_execution = Measure(
        options.iterations, options.reuse, options.retain_samples, [&]() {
        run_decode_batch();
    });

    const std::string legacy_reason = options.skip_legacy
        ? "disabled by --skip-legacy"
        : LegacyUnavailableReason(options, codec);
    const bool legacy_available = legacy_reason.empty();
    Summary legacy_encode = Summary();
    Summary legacy_decode = Summary();
    if (legacy_available)
    {
        const unsigned encode_work_count = leo_encode_work_count(options.k, options.r);
        const unsigned decode_work_count = leo_decode_work_count(options.k, options.r);
        std::vector<std::unique_ptr<LegacyStripe> > legacy;
        legacy.reserve(options.batch);
        for (size_t i = 0; i < options.batch; ++i)
        {
            legacy.push_back(std::unique_ptr<LegacyStripe>(new LegacyStripe));
            InitializeLegacyStripe(*legacy.back(), *stripes[i], options,
                encode_work_count, decode_work_count);
            RequireLegacy(leo_encode(options.bytes, options.k, options.r,
                encode_work_count, &stripes[i]->original[0], &legacy.back()->encode_work[0]),
                "legacy correctness encode");
            for (uint32_t parity = 0; parity < options.r; ++parity)
            {
                if (memcmp(legacy.back()->encode_work[parity],
                    stripes[i]->recovery_output[parity], static_cast<size_t>(options.bytes)) != 0)
                    Fail("legacy parity differs from Leopard2 legacy-high parity");
            }
            RequireLegacy(leo_decode(options.bytes, options.k, options.r,
                decode_work_count, &legacy.back()->received_original[0],
                &legacy.back()->received_recovery[0], &legacy.back()->decode_work[0]),
                "legacy correctness decode");
            for (size_t loss_i = 0; loss_i < losses.size(); ++loss_i)
            {
                const uint32_t index = losses[loss_i];
                const uint8_t* expected = stripes[i]->original_storage.bytes() +
                    static_cast<size_t>(index) * static_cast<size_t>(options.bytes);
                if (memcmp(expected, legacy.back()->decode_work[index],
                    static_cast<size_t>(options.bytes)) != 0)
                    Fail("legacy decode restored a shard incorrectly");
            }
        }
        for (size_t warmup = 0; warmup < options.warmup; ++warmup)
        {
            for (size_t i = 0; i < legacy.size(); ++i)
            {
                RequireLegacy(leo_encode(options.bytes, options.k, options.r,
                    encode_work_count, &stripes[i]->original[0], &legacy[i]->encode_work[0]),
                    "legacy encode warmup");
                RequireLegacy(leo_decode(options.bytes, options.k, options.r,
                    decode_work_count, &legacy[i]->received_original[0],
                    &legacy[i]->received_recovery[0], &legacy[i]->decode_work[0]),
                    "legacy decode warmup");
            }
        }
        legacy_encode = Measure(
            options.iterations, options.reuse, options.retain_samples, [&]() {
            for (size_t i = 0; i < legacy.size(); ++i)
                RequireLegacy(leo_encode(options.bytes, options.k, options.r,
                    encode_work_count, &stripes[i]->original[0], &legacy[i]->encode_work[0]),
                    "timed legacy encode");
        });
        legacy_decode = Measure(
            options.iterations, options.reuse, options.retain_samples, [&]() {
            for (size_t i = 0; i < legacy.size(); ++i)
                RequireLegacy(leo_decode(options.bytes, options.k, options.r,
                    decode_work_count, &legacy[i]->received_original[0],
                    &legacy[i]->received_recovery[0], &legacy[i]->decode_work[0]),
                    "timed legacy decode");
        });
    }

#if defined(LEO2_HIGH_DECODE_COPY_ATTRIBUTION)
    const HighEvaluatorAttribution high_evaluator =
        ReadHighEvaluatorAttribution(leo2_codec_field(codec), options);
#endif

    const uint64_t batch = static_cast<uint64_t>(options.batch);
    const uint64_t encode_input_bytes = CheckedU64Product(
        CheckedU64Product(options.k, options.bytes, "encode input byte count"),
        batch, "encode input byte count");
    const uint64_t encode_output_bytes = CheckedU64Product(
        CheckedU64Product(options.r, options.bytes, "encode output byte count"),
        batch, "encode output byte count");
    const uint64_t decode_input_shards = static_cast<uint64_t>(options.k - options.losses) + options.r;
    const uint64_t decode_input_bytes = CheckedU64Product(
        CheckedU64Product(decode_input_shards, options.bytes, "decode input byte count"),
        batch, "decode input byte count");
    const uint64_t decode_output_bytes = CheckedU64Product(
        CheckedU64Product(options.losses, options.bytes, "decode output byte count"),
        batch, "decode output byte count");
    const uint64_t encode_scratch_batch = CheckedU64Sum(
        CheckedU64Product(
            static_cast<uint64_t>(encode_scratch_bytes), batch,
            "batch encode scratch byte count"),
        static_cast<uint64_t>(encode_batch_preflight_bytes),
        "batch encode scratch byte count");
    const uint64_t decode_scratch_batch = CheckedU64Sum(
        CheckedU64Product(
            static_cast<uint64_t>(decode_scratch_bytes), batch,
            "batch decode scratch byte count"),
        static_cast<uint64_t>(decode_batch_preflight_bytes),
        "batch decode scratch byte count");

    std::ostringstream json;
    json << std::fixed << std::setprecision(6);
    json << "{\n"
         << "  \"schema\": \"leopard2-benchmark-v"
         << schema_version << "\",\n"
         << "  \"build\": {\n"
         << "    \"compiler\": \"" << CompilerName() << "\",\n"
         << "    \"compiler_version\": \"" << JsonEscape(CompilerVersion()) << "\",\n"
         << "    \"cplusplus\": " << __cplusplus;
#if defined(LEO2_BENCHMARK_SOURCE_ATTESTATION)
    if (options.attest_source)
    {
        json << ",\n"
             << "    \"source_commit\": \""
             << LEO2_BENCHMARK_SOURCE_COMMIT << "\",\n"
             << "    \"source_tree\": \""
             << LEO2_BENCHMARK_SOURCE_TREE << "\",\n"
             << "    \"source_tracked_dirty\": "
             << (LEO2_BENCHMARK_SOURCE_TRACKED_DIRTY ? "true" : "false");
    }
#endif
    json << "\n"
         << "  },\n"
         << "  \"parameters\": {\n"
         << "    \"K\": " << options.k << ",\n"
         << "    \"R\": " << options.r << ",\n"
         << "    \"requested_profile\": \"" << ProfileName(options.profile) << "\",\n"
         << "    \"requested_field\": \"" << FieldName(options.field) << "\",\n"
         << "    \"requested_backend\": \"" << BackendName(options.backend) << "\",\n"
         << "    \"force_generic_decode\": "
         << (options.force_generic_decode ? "true" : "false") << ",\n"
         << "    \"force_specialized_decode\": "
         << (options.force_specialized_decode ? "true" : "false") << ",\n"
         << "    \"force_tiled_decode\": "
         << (options.force_tiled_decode ? "true" : "false") << ",\n"
         << "    \"force_materialized_decode\": "
         << (options.force_materialized_decode ? "true" : "false") << ",\n";
#if defined(LEO2_HIGH_LOW_DUALITY_ATTRIBUTION)
    json << "    \"force_translated_low\": "
         << (options.force_translated_low ? "true" : "false") << ",\n"
         << "    \"force_native_high\": "
         << (options.force_native_high ? "true" : "false") << ",\n";
#endif
#if defined(LEO2_HIGH_DECODE_COPY_ATTRIBUTION)
    json << "    \"high_evaluator_mode\": \""
         << HighEvaluatorModeName(options.high_evaluator_mode) << "\",\n";
#endif
    if (extended_schema)
    {
        json << "    \"skip_legacy\": "
             << (options.skip_legacy ? "true" : "false") << ",\n"
             << "    \"retain_samples\": "
             << (options.retain_samples ? "true" : "false") << ",\n";
        if (options.report_decode_path)
            json << "    \"report_decode_path\": true,\n";
        if (options.report_direct_executor)
            json << "    \"report_direct_executor\": true,\n";
        if (options.attest_source)
            json << "    \"attest_source\": true,\n";
    }
    json << "    \"shard_bytes\": " << options.bytes << ",\n"
         << "    \"loss_count\": " << options.losses << ",\n"
         << "    \"missing_original_indices\": [";
    for (size_t i = 0; i < losses.size(); ++i)
    {
        if (i != 0)
            json << ", ";
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
         << "    \"profile\": \"" << ProfileName(leo2_codec_profile(codec)) << "\",\n"
         << "    \"field\": \"" << FieldName(leo2_codec_field(codec)) << "\",\n"
         << "    \"backend\": \"" << BackendName(leo2_context_backend(context)) << "\",\n"
         << "    \"thread_count\": " << leo2_context_thread_count(context) << ",\n"
         << "    \"parent_count\": " << leo2_codec_parent_count(codec) << ",\n"
         << "    \"padded_side\": " << leo2_codec_padded_side(codec);
    if (options.report_decode_path)
    {
        json << ",\n"
             << "    \"selected_decode_path\": \""
             << leopard2_internal::DecodePathName(decode_path_info.path)
             << "\",\n"
             << "    \"selected_decode_rule\": \""
             << leopard2_internal::DecodePathRuleName(decode_path_info.rule)
             << "\",\n"
             << "    \"decode_required_work_slots\": "
             << decode_path_info.required_work_slots << ",\n"
             << "    \"decode_aligned_prefix_bytes\": "
             << decode_path_info.aligned_prefix_bytes << ",\n"
             << "    \"decode_tail_bytes\": "
             << decode_path_info.tail_bytes << ",\n"
             << "    \"decode_rounded_bytes\": "
             << decode_path_info.rounded_shard_bytes << ",\n"
             << "    \"decode_multi_item_batch\": "
             << (decode_path_info.multi_item_batch ? "true" : "false");
        if (options.report_direct_executor)
            json << ",\n"
                 << "    \"selected_direct_executor\": \""
                 << leopard2_internal::DirectRepairExecutorName(
                        decode_path_info.direct_executor)
                 << "\"";
    }
#if defined(LEO2_HIGH_DECODE_COPY_ATTRIBUTION)
    json << ",\n"
         << "    \"high_evaluator_mode\": \""
         << HighEvaluatorModeName(options.high_evaluator_mode) << "\"";
#endif
    json << "\n"
         << "  },\n"
         << "  \"correctness\": {\n"
         << "    \"leopard2_round_trip\": true,\n"
         << "    \"legacy_comparison\": ";
    if (legacy_available)
        json << "\"matched\"\n";
    else
        json << "null\n";
    json << "  },\n";
    if (extended_schema)
    {
        json << "  \"workload_digests\": {\n"
             << "    \"algorithm\": \"fnv1a64\",\n"
             << "    \"original_data\": \"" << HexU64(original_digest) << "\",\n"
             << "    \"transmitted_parity\": \"" << HexU64(parity_digest) << "\",\n"
             << "    \"recovered_originals\": \"" << HexU64(recovered_digest) << "\"\n"
             << "  },\n";
    }
#if defined(LEO2_HIGH_DECODE_COPY_ATTRIBUTION)
    json << "  \"high_evaluator_attribution\": {\n"
         << "    \"mode\": \""
         << HighEvaluatorModeName(options.high_evaluator_mode) << "\",\n"
         << "    \"output_blocks\": " << high_evaluator.output_blocks << ",\n"
         << "    \"fft_butterfly2_out_of_place\": "
         << high_evaluator.butterfly2_out_of_place << ",\n"
         << "    \"fft_butterfly4_out_of_place\": "
         << high_evaluator.butterfly4_out_of_place << ",\n"
         << "    \"compatibility_copy_fallbacks\": "
         << high_evaluator.compatibility_copy_fallbacks << ",\n"
         << "    \"pruned_output_blocks\": "
         << high_evaluator.pruned_output_blocks << ",\n"
         << "    \"mature_output_blocks\": "
         << high_evaluator.mature_output_blocks << ",\n"
         << "    \"invariant_passed\": true\n"
         << "  },\n";
#endif
    json << "  \"memory\": {\n"
         << "    \"scratch_alignment\": " << leo2_scratch_alignment() << ",\n"
         << "    \"encode_scratch_bytes_per_stripe\": " << encode_scratch_bytes << ",\n"
         << "    \"decode_scratch_bytes_per_stripe\": " << decode_scratch_bytes << ",\n"
         << "    \"encode_scratch_bytes_batch\": " << encode_scratch_batch << ",\n"
         << "    \"decode_scratch_bytes_batch\": " << decode_scratch_batch << "\n"
         << "  },\n"
         << "  \"metrics\": {\n"
         << "    \"codec_setup\": ";
    WriteSetupSummary(json, codec_setup, 4, options.retain_samples);
    json << ",\n    \"encode_execution\": ";
    WriteSummary(json, encode_execution, encode_input_bytes, "input_GB_per_s",
        encode_output_bytes, "parity_output_GB_per_s", 4, options.retain_samples);
    json << ",\n    \"decode_plan_setup\": ";
    WriteSetupSummary(json, plan_setup, 4, options.retain_samples);
    json << ",\n    \"decode_execution\": ";
    WriteSummary(json, decode_execution, decode_input_bytes, "offered_received_GB_per_s",
        decode_output_bytes, "repaired_output_GB_per_s", 4, options.retain_samples);
    json << ",\n    \"decode_amortized_at_reuse\": ";
    WriteAmortizedDecodeSummary(json, plan_setup, decode_execution, options.reuse,
        decode_input_bytes, decode_output_bytes, 4);
    json << ",\n    \"rate_semantics\": "
         << "\"offered_received counts all non-null shard pointers supplied; "
         << "a plan may read a deterministic subset\"";
    json << "\n  },\n"
         << "  \"legacy\": {\n"
         << "    \"available\": " << (legacy_available ? "true" : "false") << ",\n"
         << "    \"unavailable_reason\": ";
    if (legacy_available)
        json << "null,\n";
    else
        json << '"' << JsonEscape(legacy_reason) << "\",\n";
    json << "    \"codec_setup\": null,\n"
         << "    \"decode_timing_includes_setup\": true,\n"
         << "    \"encode_execution\": ";
    if (legacy_available)
        WriteSummary(json, legacy_encode, encode_input_bytes, "input_GB_per_s",
            encode_output_bytes, "parity_output_GB_per_s", 4, options.retain_samples);
    else
        json << "null";
    json << ",\n    \"decode_including_setup\": ";
    if (legacy_available)
        WriteSummary(json, legacy_decode, decode_input_bytes, "offered_received_GB_per_s",
            decode_output_bytes, "repaired_output_GB_per_s", 4, options.retain_samples);
    else
        json << "null";
    json << "\n  }\n}\n";

    if (options.output == "-")
        std::cout << json.str();
    else
    {
        std::ofstream output(options.output.c_str(), std::ios::out | std::ios::trunc);
        if (!output)
            Fail("cannot open JSON output: " + options.output);
        output << json.str();
        if (!output)
            Fail("failed writing JSON output: " + options.output);
    }

    leo2_decode_plan_destroy(plan);
    leo2_codec_destroy(codec);
    leo2_context_destroy(context);
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
        std::cerr << "leopard2 benchmark: " << error.what() << std::endl;
        return 1;
    }
}
