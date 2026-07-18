/*
    Finite comparative benchmark for packed Leopard FF8 and the experimental
    native plane-sliced FF8 XOR-circuit backend.
*/

#include "../LeopardCommon.h"
#include "../LeopardFF8.h"
#include "../LeopardFF8Xor.h"
#include "../leopard.h"
#include "../leopard_ff8xor.h"

#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <new>
#include <sstream>
#include <string>
#include <vector>

namespace {

static volatile uint64_t BenchmarkSink = 0;

struct Options
{
    bool quick;
    bool csv;
    bool include_transpose;
    unsigned warmups;
    unsigned iterations;

    Options()
        : quick(false)
        , csv(false)
        , include_transpose(false)
        , warmups(2)
        , iterations(7)
    {
    }
};

struct Environment
{
    std::string compiler;
    std::string build_type;
    std::string cpu;
    std::string simd;
};

struct Timing
{
    double median_usec;
    double best_usec;

    Timing() : median_usec(0), best_usec(0) {}
};

struct Result
{
    std::string record;
    std::string backend;
    std::string operation;
    bool transpose_included;
    unsigned original_count;
    unsigned recovery_count;
    uint64_t buffer_bytes;
    unsigned loss_count;
    unsigned warmups;
    unsigned iterations;
    Timing timing;
    uint64_t input_bytes;
    uint64_t output_bytes;
    double ratio_vs_packed;
    std::string note;

    Result()
        : transpose_included(false)
        , original_count(0)
        , recovery_count(0)
        , buffer_bytes(0)
        , loss_count(0)
        , warmups(0)
        , iterations(0)
        , input_bytes(0)
        , output_bytes(0)
        , ratio_vs_packed(0)
    {
    }
};

class BufferSet
{
public:
    BufferSet() {}
    ~BufferSet() { Clear(); }

    BufferSet(const BufferSet&) = delete;
    BufferSet& operator=(const BufferSet&) = delete;

    bool Allocate(unsigned count, size_t bytes)
    {
        Clear();
        try
        {
            Buffers.reserve(count);
            for (unsigned i = 0; i < count; ++i)
            {
                uint8_t* buffer = leopard::SIMDSafeAllocate(bytes);
                if (!buffer)
                {
                    Clear();
                    return false;
                }
                Buffers.push_back(buffer);
            }
        }
        catch (const std::bad_alloc&)
        {
            Clear();
            return false;
        }
        return true;
    }

    void Clear()
    {
        for (size_t i = 0; i < Buffers.size(); ++i)
            leopard::SIMDSafeFree(Buffers[i]);
        Buffers.clear();
    }

    bool TransferFirst(unsigned count, BufferSet& destination)
    {
        if (count > Buffers.size())
            return false;
        destination.Clear();
        try
        {
            destination.Buffers.reserve(count);
        }
        catch (const std::bad_alloc&)
        {
            return false;
        }
        for (unsigned i = 0; i < count; ++i)
        {
            destination.Buffers.push_back(Buffers[i]);
            Buffers[i] = NULL;
        }
        Clear();
        return true;
    }

    uint8_t* operator[](size_t index)
    {
        return static_cast<uint8_t*>(Buffers[index]);
    }

    const uint8_t* operator[](size_t index) const
    {
        return static_cast<const uint8_t*>(Buffers[index]);
    }

    void** Data() { return Buffers.empty() ? NULL : &Buffers[0]; }
    unsigned Count() const { return static_cast<unsigned>(Buffers.size()); }

private:
    std::vector<void*> Buffers;
};

class Random
{
public:
    explicit Random(uint64_t seed) : State(seed ? seed : UINT64_C(1)) {}

    uint64_t Next()
    {
        uint64_t value = State;
        value ^= value >> 12;
        value ^= value << 25;
        value ^= value >> 27;
        State = value;
        return value * UINT64_C(2685821657736338717);
    }

private:
    uint64_t State;
};

static uint64_t TransposeWord8x8(uint64_t value)
{
    uint64_t temp = (value ^ (value >> 7)) & UINT64_C(0x00aa00aa00aa00aa);
    value ^= temp ^ (temp << 7);
    temp = (value ^ (value >> 14)) & UINT64_C(0x0000cccc0000cccc);
    value ^= temp ^ (temp << 14);
    temp = (value ^ (value >> 28)) & UINT64_C(0x00000000f0f0f0f0);
    value ^= temp ^ (temp << 28);
    return value;
}

// packed[8*g + lane] bit plane_index ==
// plane[plane_index * (bytes / 8) + g] bit lane.
static void PackedToPlane(
    const uint8_t* packed,
    uint8_t* plane,
    uint64_t bytes)
{
    const uint64_t plane_bytes = bytes / 8;
    for (uint64_t group = 0; group < plane_bytes; ++group)
    {
        uint64_t word;
        memcpy(&word, packed + group * 8, sizeof(word));
        word = TransposeWord8x8(word);
        for (unsigned bit = 0; bit < 8; ++bit)
            plane[bit * plane_bytes + group] =
                static_cast<uint8_t>(word >> (bit * 8));
    }
}

static void PlaneToPacked(
    const uint8_t* plane,
    uint8_t* packed,
    uint64_t bytes)
{
    const uint64_t plane_bytes = bytes / 8;
    for (uint64_t group = 0; group < plane_bytes; ++group)
    {
        uint64_t word = 0;
        for (unsigned bit = 0; bit < 8; ++bit)
            word |= static_cast<uint64_t>(plane[bit * plane_bytes + group])
                << (bit * 8);
        word = TransposeWord8x8(word);
        memcpy(packed + group * 8, &word, sizeof(word));
    }
}

static void FillRandom(BufferSet& buffers, uint64_t bytes, uint64_t seed)
{
    Random random(seed);
    for (unsigned shard = 0; shard < buffers.Count(); ++shard)
    {
        for (uint64_t offset = 0; offset < bytes; offset += 8)
        {
            const uint64_t value = random.Next();
            memcpy(buffers[shard] + offset, &value, sizeof(value));
        }
    }
}

static bool Equal(const void* a, const void* b, uint64_t bytes)
{
    return memcmp(a, b, static_cast<size_t>(bytes)) == 0;
}

template <typename Function>
static bool Measure(
    unsigned warmups,
    unsigned iterations,
    Function function,
    Timing& timing,
    LeopardResult& result)
{
    typedef std::chrono::steady_clock Clock;
    for (unsigned i = 0; i < warmups; ++i)
    {
        result = function();
        if (result != Leopard_Success)
            return false;
    }

    std::vector<double> samples;
    samples.reserve(iterations);
    for (unsigned i = 0; i < iterations; ++i)
    {
        const Clock::time_point begin = Clock::now();
        result = function();
        const Clock::time_point end = Clock::now();
        if (result != Leopard_Success)
            return false;
        const double usec =
            std::chrono::duration_cast<std::chrono::duration<double, std::micro> >(
                end - begin).count();
        samples.push_back(usec);
    }

    std::sort(samples.begin(), samples.end());
    timing.best_usec = samples.front();
    const size_t middle = samples.size() / 2;
    if ((samples.size() & 1) != 0)
        timing.median_usec = samples[middle];
    else
        timing.median_usec = (samples[middle - 1] + samples[middle]) * 0.5;
    return true;
}

static double Throughput(uint64_t bytes, double usec)
{
    if (usec <= 0)
        return 0;
    // Decimal MB/s: bytes/usec is numerically equal to MB/s.
    return static_cast<double>(bytes) / usec;
}

static std::string Csv(const std::string& text)
{
    std::string escaped;
    escaped.reserve(text.size() + 2);
    escaped.push_back('"');
    for (size_t i = 0; i < text.size(); ++i)
    {
        if (text[i] == '"')
            escaped.push_back('"');
        escaped.push_back(text[i]);
    }
    escaped.push_back('"');
    return escaped;
}

class Reporter
{
public:
    Reporter(bool csv, const Environment& environment)
        : CSV(csv), Env(environment)
    {
    }

    void Begin(const Options& options)
    {
        const leopard::ff8xor::CircuitStatistics multiply =
            leopard::ff8xor::GetMultiplyCircuitStatistics();
        const leopard::ff8xor::CircuitStatistics fft =
            leopard::ff8xor::GetFFTCircuitStatistics();
        const leopard::ff8xor::CircuitStatistics ifft =
            leopard::ff8xor::GetIFFTCircuitStatistics();
        if (CSV)
        {
            std::cout
                << "record,backend,operation,transpose_included,k,r,buffer_bytes,loss_count,"
                << "warmups,iterations,median_us,best_us,median_input_MBps,best_input_MBps,"
                << "median_output_MBps,best_output_MBps,ratio_vs_packed,compiler,build_type,"
                << "cpu,simd,checksum,gate_min,gate_max,gate_average,note\n";
            GateCSV("multiply",
                multiply.MinimumGateCount,
                multiply.MaximumGateCount,
                multiply.AverageGateCount);
            GateCSV("fft",
                fft.MinimumGateCount,
                fft.MaximumGateCount,
                fft.AverageGateCount);
            GateCSV("ifft",
                ifft.MinimumGateCount,
                ifft.MaximumGateCount,
                ifft.AverageGateCount);
            return;
        }

        std::cout << "FF8 XOR-circuit comparative benchmark\n"
            << "compiler: " << Env.compiler << '\n'
            << "build: " << Env.build_type << '\n'
            << "cpu: " << Env.cpu << '\n'
            << "simd: " << Env.simd << '\n'
            << "mode: " << (options.quick ? "quick" : "full")
            << ", warmups=" << options.warmups
            << ", iterations=" << options.iterations
            << ", packed-boundary transpose measurements="
            << (options.include_transpose ? "included" : "disabled") << '\n'
            << "circuit checksum: "
            << leopard::ff8xor::GetCircuitChecksum() << '\n'
            << std::fixed << std::setprecision(3)
            << "gates multiply="
            << multiply.MinimumGateCount << ".."
            << multiply.MaximumGateCount << " avg="
            << multiply.AverageGateCount
            << ", fft=" << fft.MinimumGateCount << ".."
            << fft.MaximumGateCount << " avg="
            << fft.AverageGateCount
            << ", ifft=" << ifft.MinimumGateCount << ".."
            << ifft.MaximumGateCount << " avg="
            << ifft.AverageGateCount << "\n\n";
    }

    void Print(const Result& result)
    {
        const double median_input = Throughput(
            result.input_bytes, result.timing.median_usec);
        const double best_input = Throughput(
            result.input_bytes, result.timing.best_usec);
        const double median_output = Throughput(
            result.output_bytes, result.timing.median_usec);
        const double best_output = Throughput(
            result.output_bytes, result.timing.best_usec);

        if (CSV)
        {
            std::cout << Csv(result.record) << ','
                << Csv(result.backend) << ',' << Csv(result.operation) << ','
                << (result.transpose_included ? 1 : 0) << ','
                << result.original_count << ',' << result.recovery_count << ','
                << result.buffer_bytes << ',' << result.loss_count << ','
                << result.warmups << ',' << result.iterations << ','
                << std::fixed << std::setprecision(3)
                << result.timing.median_usec << ',' << result.timing.best_usec << ','
                << median_input << ',' << best_input << ','
                << median_output << ',' << best_output << ','
                << result.ratio_vs_packed << ','
                << Csv(Env.compiler) << ',' << Csv(Env.build_type) << ','
                << Csv(Env.cpu) << ',' << Csv(Env.simd) << ','
                << Csv("") << ",,,," << Csv(result.note) << '\n';
            return;
        }

        std::cout << std::left << std::setw(25) << result.backend
            << std::setw(10) << result.operation << std::right
            << " k=" << std::setw(3) << result.original_count
            << " r=" << std::setw(3) << result.recovery_count
            << " B=" << std::setw(8) << result.buffer_bytes
            << " loss=" << std::setw(3) << result.loss_count
            << " transpose=" << (result.transpose_included ? "yes" : "no ")
            << std::fixed << std::setprecision(3)
            << " median=" << std::setw(10) << result.timing.median_usec << " us"
            << " best=" << std::setw(10) << result.timing.best_usec << " us"
            << " input=" << std::setw(10) << median_input << " MB/s"
            << " output=" << std::setw(10) << median_output << " MB/s";
        if (result.ratio_vs_packed > 0)
            std::cout << " ratio=" << result.ratio_vs_packed << 'x';
        if (!result.note.empty())
            std::cout << " (" << result.note << ')';
        std::cout << '\n';
    }

    void Skip(
        unsigned original_count,
        unsigned recovery_count,
        uint64_t buffer_bytes,
        unsigned loss_count,
        const std::string& note)
    {
        Result result;
        result.record = "skip";
        result.backend = "all";
        result.operation = "skip";
        result.original_count = original_count;
        result.recovery_count = recovery_count;
        result.buffer_bytes = buffer_bytes;
        result.loss_count = loss_count;
        result.note = note;
        Print(result);
    }

private:
    void GateCSV(const char* operation, unsigned minimum, unsigned maximum, double average)
    {
        std::cout << Csv("metadata") << ',' << Csv("generated") << ','
            << Csv(operation) << ",0,0,0,0,0,0,0,0,0,0,0,0,0,0,"
            << Csv(Env.compiler) << ',' << Csv(Env.build_type) << ','
            << Csv(Env.cpu) << ',' << Csv(Env.simd) << ','
            << Csv(leopard::ff8xor::GetCircuitChecksum()) << ','
            << minimum << ',' << maximum << ',' << std::fixed
            << std::setprecision(6) << average << ',' << Csv("") << '\n';
    }

    bool CSV;
    Environment Env;
};

static std::string CompilerName()
{
#if defined(__clang__)
    return std::string("Clang ") + __clang_version__;
#elif defined(__GNUC__)
    return std::string("GCC ") + __VERSION__;
#elif defined(_MSC_VER)
    std::ostringstream text;
    text << "MSVC " << _MSC_VER;
    return text.str();
#else
    return "unknown";
#endif
}

static std::string BuildType()
{
#if defined(NDEBUG)
    return "Release (NDEBUG)";
#elif defined(__OPTIMIZE__) || (defined(_MSC_VER) && !defined(_DEBUG))
    return "optimized (NDEBUG unset)";
#else
    return "Debug/unoptimized";
#endif
}

static std::string Trim(const std::string& input)
{
    const size_t begin = input.find_first_not_of(" \t\r\n");
    if (begin == std::string::npos)
        return "";
    const size_t end = input.find_last_not_of(" \t\r\n");
    return input.substr(begin, end - begin + 1);
}

static std::string CPUName()
{
#if defined(__linux__)
    std::ifstream cpuinfo("/proc/cpuinfo");
    std::string line;
    while (std::getline(cpuinfo, line))
    {
        if (line.compare(0, 10, "model name") == 0)
        {
            const size_t colon = line.find(':');
            if (colon != std::string::npos)
                return Trim(line.substr(colon + 1));
        }
    }
#elif defined(_WIN32)
    const char* identifier = std::getenv("PROCESSOR_IDENTIFIER");
    if (identifier)
        return identifier;
#endif
    return "unavailable";
}

static std::string SIMDName()
{
    std::string packed = "portable";
#if defined(LEO_TRY_AVX2)
    if (leopard::CpuHasAVX2)
        packed = "AVX2";
    else
#endif
    if (leopard::CpuHasSSSE3)
        packed = "128-bit/SSSE3";
#if defined(LEO_TRY_NEON)
    else if (leopard::CpuHasNeon)
        packed = "NEON";
#endif
    return std::string("packed=") + packed + "; ff8xor=" +
        leopard::ff8xor::GetKernelBackendName();
}

static Environment GetEnvironment()
{
    Environment environment;
    environment.compiler = CompilerName();
    environment.build_type = BuildType();
    environment.cpu = CPUName();
    environment.simd = SIMDName();
    return environment;
}

static Result MakeResult(
    const Options& options,
    const char* backend,
    const char* operation,
    bool transpose,
    unsigned original_count,
    unsigned recovery_count,
    uint64_t buffer_bytes,
    unsigned loss_count,
    const Timing& timing,
    uint64_t input_bytes,
    uint64_t output_bytes,
    double ratio,
    const char* note = "")
{
    Result result;
    result.record = "benchmark";
    result.backend = backend;
    result.operation = operation;
    result.transpose_included = transpose;
    result.original_count = original_count;
    result.recovery_count = recovery_count;
    result.buffer_bytes = buffer_bytes;
    result.loss_count = loss_count;
    result.warmups = options.warmups;
    result.iterations = options.iterations;
    result.timing = timing;
    result.input_bytes = input_bytes;
    result.output_bytes = output_bytes;
    result.ratio_vs_packed = ratio;
    result.note = note;
    return result;
}

static void AddUnique(std::vector<unsigned>& values, unsigned value)
{
    if (std::find(values.begin(), values.end(), value) == values.end())
        values.push_back(value);
}

static std::vector<unsigned> LossCounts(unsigned recovery_count)
{
    std::vector<unsigned> losses;
    AddUnique(losses, 1);
    AddUnique(losses, std::min(4U, recovery_count));
    AddUnique(losses, std::max(1U, recovery_count / 2));
    AddUnique(losses, recovery_count);
    return losses;
}

static std::vector<unsigned> ShuffledIndices(unsigned count, uint64_t seed)
{
    std::vector<unsigned> indices(count);
    for (unsigned i = 0; i < count; ++i)
        indices[i] = i;
    Random random(seed);
    for (unsigned i = count; i > 1; --i)
    {
        const unsigned other = static_cast<unsigned>(random.Next() % i);
        std::swap(indices[i - 1], indices[other]);
    }
    return indices;
}

static bool CheckTransposeHelper()
{
    BufferSet packed;
    BufferSet plane;
    BufferSet roundtrip;
    if (!packed.Allocate(1, 64) || !plane.Allocate(1, 64) ||
        !roundtrip.Allocate(1, 64))
        return false;
    FillRandom(packed, 64, UINT64_C(0x746573745f747270));
    PackedToPlane(packed[0], plane[0], 64);
    PlaneToPacked(plane[0], roundtrip[0], 64);
    return Equal(packed[0], roundtrip[0], 64);
}

static bool VerifyPackedRecovery(
    const BufferSet& packed_recovery,
    const BufferSet& plane_work,
    unsigned recovery_count,
    uint64_t buffer_bytes,
    uint8_t* scratch)
{
    for (unsigned i = 0; i < recovery_count; ++i)
    {
        PlaneToPacked(plane_work[i], scratch, buffer_bytes);
        if (!Equal(packed_recovery[i], scratch, buffer_bytes))
        {
            std::cerr << "encode equivalence mismatch in recovery shard "
                << i << '\n';
            return false;
        }
    }
    return true;
}

static bool RunParameter(
    const Options& options,
    Reporter& reporter,
    unsigned original_count,
    unsigned recovery_count,
    uint64_t buffer_bytes)
{
    const unsigned packed_encode_count =
        leo_encode_work_count(original_count, recovery_count);
    const unsigned xor_encode_count =
        leo_ff8xor_encode_work_count(original_count, recovery_count);
    if (packed_encode_count == 0 || xor_encode_count != packed_encode_count)
    {
        reporter.Skip(original_count, recovery_count, buffer_bytes, 0,
            "invalid or inconsistent encode work count");
        return true;
    }

    BufferSet packed_original;
    BufferSet plane_original;
    BufferSet scratch;
    if (!packed_original.Allocate(original_count, static_cast<size_t>(buffer_bytes)) ||
        !plane_original.Allocate(original_count, static_cast<size_t>(buffer_bytes)) ||
        !scratch.Allocate(1, static_cast<size_t>(buffer_bytes)))
    {
        reporter.Skip(original_count, recovery_count, buffer_bytes, 0,
            "allocation failed for originals");
        return true;
    }

    const uint64_t seed = UINT64_C(0xff8c000000000000) ^
        (static_cast<uint64_t>(original_count) << 40) ^
        (static_cast<uint64_t>(recovery_count) << 24) ^ buffer_bytes;
    FillRandom(packed_original, buffer_bytes, seed);
    for (unsigned i = 0; i < original_count; ++i)
        PackedToPlane(packed_original[i], plane_original[i], buffer_bytes);

    std::vector<const void*> packed_original_ptrs(original_count);
    std::vector<const void*> plane_original_ptrs(original_count);
    for (unsigned i = 0; i < original_count; ++i)
    {
        packed_original_ptrs[i] = packed_original[i];
        plane_original_ptrs[i] = plane_original[i];
    }

    LeopardResult api_result = Leopard_Success;
    Timing packed_encode_timing;
    BufferSet packed_encode_work;
    if (!packed_encode_work.Allocate(
            packed_encode_count, static_cast<size_t>(buffer_bytes)))
    {
        reporter.Skip(original_count, recovery_count, buffer_bytes, 0,
            "allocation failed for packed encode work");
        return true;
    }
    if (!Measure(options.warmups, options.iterations,
        [&]() -> LeopardResult {
            return leo_encode(buffer_bytes, original_count, recovery_count,
                packed_encode_count, &packed_original_ptrs[0],
                packed_encode_work.Data());
        }, packed_encode_timing, api_result))
    {
        if (api_result == Leopard_TooMuchData)
        {
            reporter.Skip(original_count, recovery_count, buffer_bytes, 0,
                "packed backend reports unsupported parameters");
            return true;
        }
        std::cerr << "packed encode failed: " << leo_result_string(api_result) << '\n';
        return false;
    }

    const uint64_t encode_input_bytes = buffer_bytes * original_count;
    const uint64_t encode_output_bytes = buffer_bytes * recovery_count;
    reporter.Print(MakeResult(options, "packed_ff8", "encode", false,
        original_count, recovery_count, buffer_bytes, 0,
        packed_encode_timing, encode_input_bytes, encode_output_bytes, 1.0));

    BufferSet packed_recovery;
    if (!packed_encode_work.TransferFirst(recovery_count, packed_recovery))
    {
        std::cerr << "unable to retain packed recovery buffers\n";
        return false;
    }

    Timing xor_encode_timing;
    BufferSet xor_encode_work;
    if (!xor_encode_work.Allocate(
            xor_encode_count, static_cast<size_t>(buffer_bytes)))
    {
        reporter.Skip(original_count, recovery_count, buffer_bytes, 0,
            "allocation failed for ff8xor encode work");
        return true;
    }
    if (!Measure(options.warmups, options.iterations,
        [&]() -> LeopardResult {
            return leo_ff8xor_encode(buffer_bytes, original_count, recovery_count,
                xor_encode_count, &plane_original_ptrs[0],
                xor_encode_work.Data());
        }, xor_encode_timing, api_result))
    {
        if (api_result == Leopard_TooMuchData)
        {
            reporter.Skip(original_count, recovery_count, buffer_bytes, 0,
                "ff8xor backend reports unsupported parameters");
            return true;
        }
        std::cerr << "ff8xor encode failed: " << leo_result_string(api_result) << '\n';
        return false;
    }
    if (!VerifyPackedRecovery(packed_recovery, xor_encode_work,
            recovery_count, buffer_bytes, scratch[0]))
        return false;

    reporter.Print(MakeResult(options, "ff8xor_native", "encode", false,
        original_count, recovery_count, buffer_bytes, 0,
        xor_encode_timing, encode_input_bytes, encode_output_bytes,
        packed_encode_timing.median_usec / xor_encode_timing.median_usec));

    if (options.include_transpose)
    {
        BufferSet packed_output;
        if (!packed_output.Allocate(
                recovery_count, static_cast<size_t>(buffer_bytes)))
        {
            reporter.Skip(original_count, recovery_count, buffer_bytes, 0,
                "allocation failed for transpose-inclusive encode output");
        }
        else
        {
            Timing included_timing;
            if (!Measure(options.warmups, options.iterations,
                [&]() -> LeopardResult {
                    for (unsigned i = 0; i < original_count; ++i)
                        PackedToPlane(packed_original[i], plane_original[i], buffer_bytes);
                    LeopardResult result = leo_ff8xor_encode(
                        buffer_bytes, original_count, recovery_count,
                        xor_encode_count, &plane_original_ptrs[0],
                        xor_encode_work.Data());
                    if (result != Leopard_Success)
                        return result;
                    for (unsigned i = 0; i < recovery_count; ++i)
                        PlaneToPacked(xor_encode_work[i], packed_output[i], buffer_bytes);
                    return Leopard_Success;
                }, included_timing, api_result))
            {
                std::cerr << "transpose-inclusive ff8xor encode failed: "
                    << leo_result_string(api_result) << '\n';
                return false;
            }
            for (unsigned i = 0; i < recovery_count; ++i)
            {
                if (!Equal(packed_recovery[i], packed_output[i], buffer_bytes))
                {
                    std::cerr << "transpose-inclusive encode mismatch at shard "
                        << i << '\n';
                    return false;
                }
            }
            reporter.Print(MakeResult(options, "ff8xor_packed_boundary", "encode", true,
                original_count, recovery_count, buffer_bytes, 0,
                included_timing, encode_input_bytes, encode_output_bytes,
                packed_encode_timing.median_usec / included_timing.median_usec));
        }
    }

    BufferSet plane_recovery;
    if (!xor_encode_work.TransferFirst(recovery_count, plane_recovery))
    {
        std::cerr << "unable to retain ff8xor recovery buffers\n";
        return false;
    }

    const std::vector<unsigned> loss_counts = LossCounts(recovery_count);
    for (size_t loss_case = 0; loss_case < loss_counts.size(); ++loss_case)
    {
        const unsigned loss_count = loss_counts[loss_case];
        const uint64_t loss_seed = seed ^
            (static_cast<uint64_t>(loss_count) << 8) ^ UINT64_C(0xd3c0de);
        const std::vector<unsigned> original_order =
            ShuffledIndices(original_count, loss_seed);
        const std::vector<unsigned> recovery_order =
            ShuffledIndices(recovery_count, loss_seed ^ UINT64_C(0x7265636f76657279));

        std::vector<bool> original_missing(original_count, false);
        std::vector<bool> recovery_available(recovery_count, false);
        for (unsigned i = 0; i < loss_count; ++i)
        {
            original_missing[original_order[i]] = true;
            recovery_available[recovery_order[i]] = true;
        }

        std::vector<const void*> packed_decode_original(original_count);
        std::vector<const void*> plane_decode_original(original_count);
        std::vector<const void*> packed_decode_recovery(recovery_count);
        std::vector<const void*> plane_decode_recovery(recovery_count);
        for (unsigned i = 0; i < original_count; ++i)
        {
            packed_decode_original[i] = original_missing[i] ? NULL : packed_original[i];
            plane_decode_original[i] = original_missing[i] ? NULL : plane_original[i];
        }
        for (unsigned i = 0; i < recovery_count; ++i)
        {
            packed_decode_recovery[i] = recovery_available[i]
                ? packed_recovery[i] : NULL;
            plane_decode_recovery[i] = recovery_available[i]
                ? plane_recovery[i] : NULL;
        }

        const unsigned packed_decode_count =
            leo_decode_work_count(original_count, recovery_count);
        const unsigned xor_decode_count =
            leo_ff8xor_decode_work_count(original_count, recovery_count);
        if (packed_decode_count == 0 || xor_decode_count != packed_decode_count)
        {
            reporter.Skip(original_count, recovery_count, buffer_bytes, loss_count,
                "invalid or inconsistent decode work count");
            continue;
        }

        BufferSet packed_decode_work;
        if (!packed_decode_work.Allocate(
                packed_decode_count, static_cast<size_t>(buffer_bytes)))
        {
            reporter.Skip(original_count, recovery_count, buffer_bytes, loss_count,
                "allocation failed for packed decode work");
            continue;
        }
        Timing packed_decode_timing;
        if (!Measure(options.warmups, options.iterations,
            [&]() -> LeopardResult {
                return leo_decode(buffer_bytes, original_count, recovery_count,
                    packed_decode_count, &packed_decode_original[0],
                    &packed_decode_recovery[0], packed_decode_work.Data());
            }, packed_decode_timing, api_result))
        {
            std::cerr << "packed decode failed: " << leo_result_string(api_result) << '\n';
            return false;
        }
        for (unsigned i = 0; i < loss_count; ++i)
        {
            const unsigned index = original_order[i];
            if (!Equal(packed_original[index], packed_decode_work[index], buffer_bytes))
            {
                std::cerr << "packed decode mismatch at original shard " << index << '\n';
                return false;
            }
        }

        const uint64_t decode_input_bytes = buffer_bytes * original_count;
        const uint64_t decode_output_bytes = buffer_bytes * loss_count;
        reporter.Print(MakeResult(options, "packed_ff8", "decode", false,
            original_count, recovery_count, buffer_bytes, loss_count,
            packed_decode_timing, decode_input_bytes, decode_output_bytes, 1.0));
        packed_decode_work.Clear();

        BufferSet xor_decode_work;
        if (!xor_decode_work.Allocate(
                xor_decode_count, static_cast<size_t>(buffer_bytes)))
        {
            reporter.Skip(original_count, recovery_count, buffer_bytes, loss_count,
                "allocation failed for ff8xor decode work");
            continue;
        }
        Timing xor_decode_timing;
        if (!Measure(options.warmups, options.iterations,
            [&]() -> LeopardResult {
                return leo_ff8xor_decode(buffer_bytes, original_count, recovery_count,
                    xor_decode_count, &plane_decode_original[0],
                    &plane_decode_recovery[0], xor_decode_work.Data());
            }, xor_decode_timing, api_result))
        {
            std::cerr << "ff8xor decode failed: " << leo_result_string(api_result) << '\n';
            return false;
        }
        for (unsigned i = 0; i < loss_count; ++i)
        {
            const unsigned index = original_order[i];
            if (!Equal(plane_original[index], xor_decode_work[index], buffer_bytes))
            {
                std::cerr << "native ff8xor decode mismatch at original shard "
                    << index << '\n';
                return false;
            }
            PlaneToPacked(xor_decode_work[index], scratch[0], buffer_bytes);
            if (!Equal(packed_original[index], scratch[0], buffer_bytes))
            {
                std::cerr << "packed ff8xor decode equivalence mismatch at shard "
                    << index << '\n';
                return false;
            }
        }
        reporter.Print(MakeResult(options, "ff8xor_native", "decode", false,
            original_count, recovery_count, buffer_bytes, loss_count,
            xor_decode_timing, decode_input_bytes, decode_output_bytes,
            packed_decode_timing.median_usec / xor_decode_timing.median_usec));

        if (options.include_transpose)
        {
            BufferSet packed_output;
            if (!packed_output.Allocate(
                    loss_count, static_cast<size_t>(buffer_bytes)))
            {
                reporter.Skip(original_count, recovery_count, buffer_bytes, loss_count,
                    "allocation failed for transpose-inclusive decode output");
            }
            else
            {
                Timing included_timing;
                if (!Measure(options.warmups, options.iterations,
                    [&]() -> LeopardResult {
                        for (unsigned i = 0; i < original_count; ++i)
                            if (!original_missing[i])
                                PackedToPlane(packed_original[i], plane_original[i],
                                    buffer_bytes);
                        for (unsigned i = 0; i < recovery_count; ++i)
                            if (recovery_available[i])
                                PackedToPlane(packed_recovery[i], plane_recovery[i],
                                    buffer_bytes);
                        LeopardResult result = leo_ff8xor_decode(
                            buffer_bytes, original_count, recovery_count,
                            xor_decode_count, &plane_decode_original[0],
                            &plane_decode_recovery[0], xor_decode_work.Data());
                        if (result != Leopard_Success)
                            return result;
                        for (unsigned i = 0; i < loss_count; ++i)
                        {
                            const unsigned index = original_order[i];
                            PlaneToPacked(xor_decode_work[index], packed_output[i],
                                buffer_bytes);
                        }
                        return Leopard_Success;
                    }, included_timing, api_result))
                {
                    std::cerr << "transpose-inclusive ff8xor decode failed: "
                        << leo_result_string(api_result) << '\n';
                    return false;
                }
                for (unsigned i = 0; i < loss_count; ++i)
                {
                    const unsigned index = original_order[i];
                    if (!Equal(packed_original[index], packed_output[i], buffer_bytes))
                    {
                        std::cerr << "transpose-inclusive decode mismatch at shard "
                            << index << '\n';
                        return false;
                    }
                }
                reporter.Print(MakeResult(options, "ff8xor_packed_boundary", "decode", true,
                    original_count, recovery_count, buffer_bytes, loss_count,
                    included_timing, decode_input_bytes, decode_output_bytes,
                    packed_decode_timing.median_usec / included_timing.median_usec));
            }
        }
    }

    if (!options.csv)
        std::cout << '\n';
    return true;
}

static void PrintMicro(
    const Options& options,
    Reporter& reporter,
    const char* backend,
    const char* operation,
    const char* note,
    const Timing& timing,
    uint64_t buffer_bytes,
    uint64_t input_bytes,
    uint64_t output_bytes)
{
    Result result = MakeResult(options, backend, operation, false,
        0, 0, buffer_bytes, 0, timing, input_bytes, output_bytes, 0, note);
    result.record = "microbenchmark";
    reporter.Print(result);
}

static bool RunMicrobenchmarks(const Options& base_options, Reporter& reporter)
{
    Options options = base_options;
    options.warmups = base_options.quick ? 2 : 3;
    options.iterations = base_options.quick ? 15 : 31;
    const uint64_t bytes = base_options.quick ? 64 * 1024 : 1024 * 1024;

    BufferSet a;
    BufferSet b;
    BufferSet c;
    if (!a.Allocate(1, static_cast<size_t>(bytes)) ||
        !b.Allocate(1, static_cast<size_t>(bytes)) ||
        !c.Allocate(1, static_cast<size_t>(bytes)))
    {
        reporter.Skip(0, 0, bytes, 0, "allocation failed for microbenchmarks");
        return true;
    }
    FillRandom(a, bytes, UINT64_C(0x6d6963726f5f6131));
    FillRandom(b, bytes, UINT64_C(0x6d6963726f5f6232));

    LeopardResult result = Leopard_Success;
    Timing timing;
    if (!Measure(options.warmups, options.iterations,
        [&]() -> LeopardResult {
            leopard::ff8::ExperimentalPackedMulAdd(bytes, a[0], b[0], 1);
            BenchmarkSink ^= a[0][0];
            return Leopard_Success;
        }, timing, result))
        return false;
    PrintMicro(options, reporter, "packed_ff8", "multiply_add",
        "log=1; existing packed multiply-add", timing, bytes, bytes * 2, bytes);

    timing = Timing();
    if (!Measure(options.warmups, options.iterations,
        [&]() -> LeopardResult {
            leopard::ff8xor::MultiplyBuffer(bytes, c[0], b[0], 1);
            BenchmarkSink ^= c[0][0];
            return Leopard_Success;
        }, timing, result))
        return false;
    PrintMicro(options, reporter, "ff8xor_native", "multiply",
        "log=1; generated whole-buffer kernel", timing, bytes, bytes, bytes);

    timing = Timing();
    if (!Measure(options.warmups, options.iterations,
        [&]() -> LeopardResult {
            leopard::ff8xor::FFTButterflyBuffer(bytes, a[0], b[0], 1);
            BenchmarkSink ^= a[0][0];
            return Leopard_Success;
        }, timing, result))
        return false;
    PrintMicro(options, reporter, "ff8xor_native", "fft_butterfly",
        "skew=1; generated two-buffer kernel", timing, bytes,
        bytes * 2, bytes * 2);

    timing = Timing();
    if (!Measure(options.warmups, options.iterations,
        [&]() -> LeopardResult {
            leopard::ff8xor::IFFTButterflyBuffer(bytes, a[0], b[0], 1);
            BenchmarkSink ^= b[0][0];
            return Leopard_Success;
        }, timing, result))
        return false;
    PrintMicro(options, reporter, "ff8xor_native", "ifft_butterfly",
        "skew=1; generated two-buffer kernel", timing, bytes,
        bytes * 2, bytes * 2);
    timing = Timing();
    if (!Measure(options.warmups, options.iterations,
        [&]() -> LeopardResult {
            PackedToPlane(a[0], c[0], bytes);
            BenchmarkSink ^= c[0][0];
            return Leopard_Success;
        }, timing, result))
        return false;
    PrintMicro(options, reporter, "transpose", "packed_to_plane",
        "8x8 bit transpose", timing, bytes, bytes, bytes);

    timing = Timing();
    if (!Measure(options.warmups, options.iterations,
        [&]() -> LeopardResult {
            PlaneToPacked(c[0], a[0], bytes);
            BenchmarkSink ^= a[0][0];
            return Leopard_Success;
        }, timing, result))
        return false;
    PrintMicro(options, reporter, "transpose", "plane_to_packed",
        "8x8 bit transpose", timing, bytes, bytes, bytes);

    if (!base_options.csv)
        std::cout << '\n';
    return true;
}

static void Usage(const char* program)
{
    std::cout << "Usage: " << program
        << " [--quick] [--csv] [--include-transpose]\n"
        << "  --quick              Run a bounded development/CI subset.\n"
        << "  --csv                Emit machine-readable CSV.\n"
        << "  --include-transpose  Also time packed boundary transposes around ff8xor.\n";
}

static bool ParseOptions(int argc, char** argv, Options& options)
{
    for (int i = 1; i < argc; ++i)
    {
        const std::string argument(argv[i]);
        if (argument == "--quick")
            options.quick = true;
        else if (argument == "--csv")
            options.csv = true;
        else if (argument == "--include-transpose")
            options.include_transpose = true;
        else if (argument == "--help" || argument == "-h")
        {
            Usage(argv[0]);
            std::exit(0);
        }
        else
        {
            std::cerr << "Unknown option: " << argument << '\n';
            Usage(argv[0]);
            return false;
        }
    }
    if (options.quick)
    {
        options.warmups = 1;
        options.iterations = 3;
    }
    return true;
}

struct Counts
{
    unsigned original;
    unsigned recovery;
};

} // namespace

int main(int argc, char** argv)
{
    Options options;
    if (!ParseOptions(argc, argv, options))
        return 2;

    const LeopardResult init_result = static_cast<LeopardResult>(leo_init());
    if (init_result != Leopard_Success)
    {
        std::cerr << "leo_init failed: " << leo_result_string(init_result) << '\n';
        return 1;
    }

    Reporter reporter(options.csv, GetEnvironment());
    reporter.Begin(options);

    if (!CheckTransposeHelper())
    {
        std::cerr << "8x8 transpose self-test failed\n";
        return 1;
    }

    if (!RunMicrobenchmarks(options, reporter))
        return 1;

    static const Counts full_counts[] = {
        { 8, 2 }, { 16, 4 }, { 32, 8 },
        { 64, 16 }, { 128, 32 }, { 128, 128 }
    };
    static const Counts quick_counts[] = {
        { 8, 2 }, { 32, 8 }
    };
    static const uint64_t full_bytes[] = {
        1024, 4096, 64 * 1024, 1024 * 1024
    };
    static const uint64_t quick_bytes[] = {
        1024, 64 * 1024
    };

    const Counts* counts = options.quick ? quick_counts : full_counts;
    const size_t counts_size = options.quick
        ? sizeof(quick_counts) / sizeof(quick_counts[0])
        : sizeof(full_counts) / sizeof(full_counts[0]);
    const uint64_t* sizes = options.quick ? quick_bytes : full_bytes;
    const size_t sizes_count = options.quick
        ? sizeof(quick_bytes) / sizeof(quick_bytes[0])
        : sizeof(full_bytes) / sizeof(full_bytes[0]);

    try
    {
        for (size_t pair = 0; pair < counts_size; ++pair)
        {
            for (size_t size = 0; size < sizes_count; ++size)
            {
                if (!RunParameter(options, reporter,
                        counts[pair].original, counts[pair].recovery,
                        sizes[size]))
                    return 1;
            }
        }
    }
    catch (const std::bad_alloc&)
    {
        std::cerr << "allocation failed while constructing benchmark metadata\n";
        return 1;
    }

    if (!options.csv)
        std::cout << "Benchmark complete. sink=" << BenchmarkSink << '\n';
    return 0;
}
