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

#include "benchmark_ff8xor_xor3.h"

#if defined(LEO_FF8XOR_XOR3_AVX512)

#include <algorithm>
#include <chrono>
#include <immintrin.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vector>

#if defined(_WIN32)
# include <malloc.h>
#endif


#if defined(_MSC_VER)
# define LEO_XOR3_NOINLINE __declspec(noinline)
# define LEO_XOR3_USED
#else
# define LEO_XOR3_NOINLINE __attribute__((noinline))
# define LEO_XOR3_USED __attribute__((used))
#endif

namespace {


volatile uint64_t BenchmarkSink = 0;

static inline uint64_t ReadTsc()
{
    // LFENCE + RDTSC is available throughout x86-64 and does not rely on the
    // independently enumerated RDTSCP feature.  Fence both sides so neither
    // the measured kernel nor surrounding bookkeeping overlaps the sample.
    _mm_lfence();
    const uint64_t value = __rdtsc();
    _mm_lfence();
    return value;
}

static inline uint64_t ReadRdpru(unsigned selector)
{
#if defined(__x86_64__) || defined(__i386__)
    uint32_t eax, edx;
    __asm__ volatile(
        "lfence\n\t"
        ".byte 0x0f, 0x01, 0xfd\n\t"
        "lfence"
        : "=a"(eax), "=d"(edx)
        : "c"(selector)
        : "memory");
    return (static_cast<uint64_t>(edx) << 32) | eax;
#else
    (void)selector;
    return 0;
#endif
}

static uint64_t Mix64(uint64_t& state)
{
    state += 0x9e3779b97f4a7c15ULL;
    uint64_t x = state;
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
    x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
    return x ^ (x >> 31);
}

#if defined(__GNUC__)
# pragma GCC diagnostic push
# pragma GCC diagnostic ignored "-Wignored-attributes"
#endif

template <typename Value>
struct VectorTraits;

template <>
struct VectorTraits<__m128i>
{
    static const unsigned kBytes = 16;
    static const unsigned kLanes = 4;

    static inline __m128i Load(const void* pointer)
    {
        return _mm_loadu_si128(reinterpret_cast<const __m128i*>(pointer));
    }

    static inline void Store(void* pointer, __m128i value)
    {
        _mm_storeu_si128(reinterpret_cast<__m128i*>(pointer), value);
    }

    static inline __m128i Set1(uint32_t value)
    {
        return _mm_set1_epi32(static_cast<int>(value));
    }
};

template <>
struct VectorTraits<__m256i>
{
    static const unsigned kBytes = 32;
    static const unsigned kLanes = 8;

    static inline __m256i Load(const void* pointer)
    {
        return _mm256_loadu_si256(reinterpret_cast<const __m256i*>(pointer));
    }

    static inline void Store(void* pointer, __m256i value)
    {
        _mm256_storeu_si256(reinterpret_cast<__m256i*>(pointer), value);
    }

    static inline __m256i Set1(uint32_t value)
    {
        return _mm256_set1_epi32(static_cast<int>(value));
    }
};

template <>
struct VectorTraits<__m512i>
{
    static const unsigned kBytes = 64;
    static const unsigned kLanes = 16;

    static inline __m512i Load(const void* pointer)
    {
        return _mm512_loadu_si512(pointer);
    }

    static inline void Store(void* pointer, __m512i value)
    {
        _mm512_storeu_si512(pointer, value);
    }

    static inline __m512i Set1(uint32_t value)
    {
        return _mm512_set1_epi32(static_cast<int>(value));
    }
};

#if defined(__GNUC__)
# pragma GCC diagnostic pop
#endif

struct ForcedXor2
{
    static const char* Name() { return "forced-xor2"; }

    static inline __m128i Apply(__m128i a, __m128i b, __m128i c)
    {
        __asm__ volatile(
            "vpxor %[b], %[a], %[a]\n\t"
            "vpxor %[c], %[a], %[a]"
            : [a] "+&x"(a)
            : [b] "x"(b), [c] "x"(c));
        return a;
    }

    static inline __m256i Apply(__m256i a, __m256i b, __m256i c)
    {
        // Inline assembly is intentional: GCC otherwise combines the two
        // intrinsics below into VPTERNLOGD, which would not be an XOR2 control.
        __asm__ volatile(
            "vpxor %[b], %[a], %[a]\n\t"
            "vpxor %[c], %[a], %[a]"
            : [a] "+&x"(a)
            : [b] "x"(b), [c] "x"(c));
        return a;
    }

    static inline __m512i Apply(__m512i a, __m512i b, __m512i c)
    {
        __asm__ volatile(
            "vpxord %[b], %[a], %[a]\n\t"
            "vpxord %[c], %[a], %[a]"
            : [a] "+&v"(a)
            : [b] "v"(b), [c] "v"(c));
        return a;
    }
};

struct AutoXor2
{
    static const char* Name() { return "auto-xor2-source"; }

    static inline __m128i Apply(__m128i a, __m128i b, __m128i c)
    {
        return _mm_xor_si128(_mm_xor_si128(a, b), c);
    }

    static inline __m256i Apply(__m256i a, __m256i b, __m256i c)
    {
        return _mm256_xor_si256(_mm256_xor_si256(a, b), c);
    }

    static inline __m512i Apply(__m512i a, __m512i b, __m512i c)
    {
        return _mm512_xor_si512(_mm512_xor_si512(a, b), c);
    }
};

struct ExplicitXor3
{
    static const char* Name() { return "explicit-vpternlog"; }

    static inline __m128i Apply(__m128i a, __m128i b, __m128i c)
    {
        return _mm_ternarylogic_epi32(a, b, c, 0x96);
    }

    static inline __m256i Apply(__m256i a, __m256i b, __m256i c)
    {
        return _mm256_ternarylogic_epi32(a, b, c, 0x96);
    }

    static inline __m512i Apply(__m512i a, __m512i b, __m512i c)
    {
        return _mm512_ternarylogic_epi32(a, b, c, 0x96);
    }
};

} // namespace


// Stable, noinline symbols make the intended lowering and code size directly
// inspectable.  The benchmark uses the inline forms above; correctness also
// exercises these probes so link-time dead stripping cannot remove them.
extern "C" LEO_XOR3_NOINLINE LEO_XOR3_USED __m128i
ff8xor_xor3_probe_forced_xmm(__m128i a, __m128i b, __m128i c)
{
    __asm__ volatile(
        "vpxor %[b], %[a], %[a]\n\t"
        "vpxor %[c], %[a], %[a]"
        : [a] "+&x"(a)
        : [b] "x"(b), [c] "x"(c));
    return a;
}

extern "C" LEO_XOR3_NOINLINE LEO_XOR3_USED __m256i
ff8xor_xor3_probe_forced_ymm(__m256i a, __m256i b, __m256i c)
{
    __asm__ volatile(
        "vpxor %[b], %[a], %[a]\n\t"
        "vpxor %[c], %[a], %[a]"
        : [a] "+&x"(a)
        : [b] "x"(b), [c] "x"(c));
    return a;
}

extern "C" LEO_XOR3_NOINLINE LEO_XOR3_USED __m512i
ff8xor_xor3_probe_forced_zmm(__m512i a, __m512i b, __m512i c)
{
    __asm__ volatile(
        "vpxord %[b], %[a], %[a]\n\t"
        "vpxord %[c], %[a], %[a]"
        : [a] "+&v"(a)
        : [b] "v"(b), [c] "v"(c));
    return a;
}

extern "C" LEO_XOR3_NOINLINE LEO_XOR3_USED __m128i
ff8xor_xor3_probe_explicit_xmm(__m128i a, __m128i b, __m128i c)
{
    return _mm_ternarylogic_epi32(a, b, c, 0x96);
}

extern "C" LEO_XOR3_NOINLINE LEO_XOR3_USED __m256i
ff8xor_xor3_probe_explicit_ymm(__m256i a, __m256i b, __m256i c)
{
    return _mm256_ternarylogic_epi32(a, b, c, 0x96);
}

extern "C" LEO_XOR3_NOINLINE LEO_XOR3_USED __m512i
ff8xor_xor3_probe_explicit_zmm(__m512i a, __m512i b, __m512i c)
{
    return _mm512_ternarylogic_epi32(a, b, c, 0x96);
}


namespace {


struct ProbeForcedXor2
{
    static const char* Name() { return "probe-forced-xor2"; }
    static inline __m128i Apply(__m128i a, __m128i b, __m128i c)
    {
        return ff8xor_xor3_probe_forced_xmm(a, b, c);
    }
    static inline __m256i Apply(__m256i a, __m256i b, __m256i c)
    {
        return ff8xor_xor3_probe_forced_ymm(a, b, c);
    }
    static inline __m512i Apply(__m512i a, __m512i b, __m512i c)
    {
        return ff8xor_xor3_probe_forced_zmm(a, b, c);
    }
};

struct ProbeExplicitXor3
{
    static const char* Name() { return "probe-explicit-vpternlog"; }
    static inline __m128i Apply(__m128i a, __m128i b, __m128i c)
    {
        return ff8xor_xor3_probe_explicit_xmm(a, b, c);
    }
    static inline __m256i Apply(__m256i a, __m256i b, __m256i c)
    {
        return ff8xor_xor3_probe_explicit_ymm(a, b, c);
    }
    static inline __m512i Apply(__m512i a, __m512i b, __m512i c)
    {
        return ff8xor_xor3_probe_explicit_zmm(a, b, c);
    }
};


static inline void OpaqueOperands(__m128i& a, __m128i& b, __m128i& c)
{
    __asm__ volatile("" : "+x"(a), "+x"(b), "+x"(c));
}

static inline void OpaqueOperands(__m256i& a, __m256i& b, __m256i& c)
{
    // Tell the optimizer that each loop-carried result is opaque without
    // emitting an instruction.  Otherwise it can replace an even number of
    // constant XOR iterations with no work, which is not a latency test.
    __asm__ volatile("" : "+x"(a), "+x"(b), "+x"(c));
}

static inline void OpaqueOperands(__m512i& a, __m512i& b, __m512i& c)
{
    __asm__ volatile("" : "+v"(a), "+v"(b), "+v"(c));
}

template <typename Value, typename Operation>
LEO_XOR3_NOINLINE uint64_t LatencyKernel(uint64_t iterations, const void*)
{
    Value a = VectorTraits<Value>::Set1(0x12345678U);
    Value b = VectorTraits<Value>::Set1(0xa5a5f00dU);
    Value c = VectorTraits<Value>::Set1(0x5a5a33ccU);
    for (uint64_t i = 0; i < iterations; ++i)
    {
        a = Operation::Apply(a, b, c);
        OpaqueOperands(a, b, c);
    }

    alignas(64) uint64_t words[8];
    VectorTraits<Value>::Store(words, a);
    return words[0] ^ words[VectorTraits<Value>::kBytes / 8 - 1];
}

template <typename Value, typename Operation>
LEO_XOR3_NOINLINE uint64_t Parallel4Kernel(uint64_t iterations, const void*)
{
    Value a0 = VectorTraits<Value>::Set1(0x01234567U);
    Value a1 = VectorTraits<Value>::Set1(0x89abcdefU);
    Value a2 = VectorTraits<Value>::Set1(0xfedcba98U);
    Value a3 = VectorTraits<Value>::Set1(0x76543210U);
    Value b0 = VectorTraits<Value>::Set1(0x11112222U);
    Value b1 = VectorTraits<Value>::Set1(0x33334444U);
    Value b2 = VectorTraits<Value>::Set1(0x55556666U);
    Value b3 = VectorTraits<Value>::Set1(0x77778888U);
    Value c0 = VectorTraits<Value>::Set1(0x9999aaaaU);
    Value c1 = VectorTraits<Value>::Set1(0xbbbbccccU);
    Value c2 = VectorTraits<Value>::Set1(0xddddeeeeU);
    Value c3 = VectorTraits<Value>::Set1(0xffff0000U);
    for (uint64_t i = 0; i < iterations; ++i)
    {
        a0 = Operation::Apply(a0, b0, c0);
        OpaqueOperands(a0, b0, c0);
        a1 = Operation::Apply(a1, b1, c1);
        OpaqueOperands(a1, b1, c1);
        a2 = Operation::Apply(a2, b2, c2);
        OpaqueOperands(a2, b2, c2);
        a3 = Operation::Apply(a3, b3, c3);
        OpaqueOperands(a3, b3, c3);
    }

    alignas(64) uint64_t words[8];
    VectorTraits<Value>::Store(words, a0);
    uint64_t result = words[0];
    VectorTraits<Value>::Store(words, a1);
    result ^= words[0];
    VectorTraits<Value>::Store(words, a2);
    result ^= words[0];
    VectorTraits<Value>::Store(words, a3);
    return result ^ words[0];
}

struct StreamContext
{
    uint8_t* a;
    uint8_t* b;
    uint8_t* c;
    uint8_t* output;
    uint64_t array_bytes;
};

template <typename Value, typename Operation>
LEO_XOR3_NOINLINE uint64_t StreamKernel(
    uint64_t sweeps,
    const void* context_void)
{
    const StreamContext& context =
        *reinterpret_cast<const StreamContext*>(context_void);
    for (uint64_t sweep = 0; sweep < sweeps; ++sweep)
    {
        for (uint64_t offset = 0;
            offset < context.array_bytes;
            offset += VectorTraits<Value>::kBytes)
        {
            Value a = VectorTraits<Value>::Load(context.a + offset);
            const Value b = VectorTraits<Value>::Load(context.b + offset);
            const Value c = VectorTraits<Value>::Load(context.c + offset);
            a = Operation::Apply(a, b, c);
            VectorTraits<Value>::Store(context.output + offset, a);
        }
    }

    uint64_t result = 0;
    memcpy(&result, context.output, sizeof(result));
    return result;
}

typedef uint64_t (*KernelFunction)(uint64_t, const void*);

struct TimedSample
{
    double nanoseconds;
    double tsc_ghz;
    double aperf_ghz;
};

static TimedSample TimeKernel(
    KernelFunction function,
    uint64_t iterations,
    const void* context,
    bool have_rdpru)
{
    const uint64_t aperf0 = have_rdpru ? ReadRdpru(1) : 0;
    const uint64_t tsc0 = ReadTsc();
    const std::chrono::steady_clock::time_point start =
        std::chrono::steady_clock::now();
    BenchmarkSink ^= function(iterations, context);
    const std::chrono::steady_clock::time_point stop =
        std::chrono::steady_clock::now();
    const uint64_t tsc1 = ReadTsc();
    const uint64_t aperf1 = have_rdpru ? ReadRdpru(1) : 0;

    TimedSample sample;
    sample.nanoseconds = static_cast<double>(
        std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start).count());
    sample.tsc_ghz = sample.nanoseconds > 0 ?
        static_cast<double>(tsc1 - tsc0) / sample.nanoseconds : 0;
    sample.aperf_ghz = have_rdpru && sample.nanoseconds > 0 ?
        static_cast<double>(aperf1 - aperf0) / sample.nanoseconds : 0;
    return sample;
}

static double Median(std::vector<double> values)
{
    std::sort(values.begin(), values.end());
    const size_t middle = values.size() / 2;
    if ((values.size() & 1) != 0)
        return values[middle];
    return (values[middle - 1] + values[middle]) * 0.5;
}

struct Summary
{
    double median_ns;
    double best_ns;
    double median_tsc_ghz;
    double median_aperf_ghz;
};

static Summary Summarize(const std::vector<TimedSample>& samples)
{
    std::vector<double> time;
    std::vector<double> tsc;
    std::vector<double> aperf;
    double best = samples[0].nanoseconds;
    for (size_t i = 0; i < samples.size(); ++i)
    {
        time.push_back(samples[i].nanoseconds);
        tsc.push_back(samples[i].tsc_ghz);
        if (samples[i].aperf_ghz > 0)
            aperf.push_back(samples[i].aperf_ghz);
        if (samples[i].nanoseconds < best)
            best = samples[i].nanoseconds;
    }

    Summary summary;
    summary.median_ns = Median(time);
    summary.best_ns = best;
    summary.median_tsc_ghz = Median(tsc);
    summary.median_aperf_ghz = aperf.empty() ? 0 : Median(aperf);
    return summary;
}

static void MeasurePair(
    const char* width,
    const char* shape,
    uint64_t working_set,
    KernelFunction baseline,
    const char* baseline_name,
    KernelFunction candidate,
    const char* candidate_name,
    uint64_t iterations,
    double operations_per_iteration,
    double logical_bytes_per_iteration,
    unsigned rounds,
    const void* context,
    bool have_rdpru)
{
    // Warm both paths outside the timed samples.
    BenchmarkSink ^= baseline(iterations, context);
    BenchmarkSink ^= candidate(iterations, context);

    std::vector<TimedSample> baseline_samples;
    std::vector<TimedSample> candidate_samples;
    for (unsigned round = 0; round < rounds; ++round)
    {
        // A-B-B-A makes each implementation see both positions in the pair.
        baseline_samples.push_back(TimeKernel(
            baseline, iterations, context, have_rdpru));
        candidate_samples.push_back(TimeKernel(
            candidate, iterations, context, have_rdpru));
        candidate_samples.push_back(TimeKernel(
            candidate, iterations, context, have_rdpru));
        baseline_samples.push_back(TimeKernel(
            baseline, iterations, context, have_rdpru));
    }

    const Summary a = Summarize(baseline_samples);
    const Summary b = Summarize(candidate_samples);
    const double operations = operations_per_iteration * iterations;
    const double logical_bytes = logical_bytes_per_iteration * iterations;
    const double a_gops = operations / a.median_ns;
    const double b_gops = operations / b.median_ns;
    const double a_gib = logical_bytes > 0 ?
        logical_bytes / a.median_ns * 1.e9 / (1024. * 1024. * 1024.) : 0;
    const double b_gib = logical_bytes > 0 ?
        logical_bytes / b.median_ns * 1.e9 / (1024. * 1024. * 1024.) : 0;

    printf("%-4s %-10s ws=%-9llu %-20s median=%10.0f ns best=%10.0f "
           "Gxor3/s=%7.2f logical_GiB/s=%7.2f traffic_GiB/s=%7.2f "
           "TSC_GHz=%5.3f APERF_GHz=%5.3f\n",
        width, shape, static_cast<unsigned long long>(working_set), baseline_name,
        a.median_ns, a.best_ns, a_gops, a_gib, a_gib * 4.,
        a.median_tsc_ghz, a.median_aperf_ghz);
    printf("%-4s %-10s ws=%-9llu %-20s median=%10.0f ns best=%10.0f "
           "Gxor3/s=%7.2f logical_GiB/s=%7.2f traffic_GiB/s=%7.2f "
           "TSC_GHz=%5.3f APERF_GHz=%5.3f "
           "speedup=%6.3fx\n",
        width, shape, static_cast<unsigned long long>(working_set), candidate_name,
        b.median_ns, b.best_ns, b_gops, b_gib, b_gib * 4.,
        b.median_tsc_ghz, b.median_aperf_ghz,
        a.median_ns / b.median_ns);
}

template <typename Value, typename Operation>
static bool CheckOperation(const char* width)
{
    alignas(64) uint32_t input_a[16];
    alignas(64) uint32_t input_b[16];
    alignas(64) uint32_t input_c[16];
    alignas(64) uint32_t output[16];

    for (unsigned truth = 0; truth < 8; ++truth)
    {
        const uint32_t a = (truth & 1) ? ~0U : 0;
        const uint32_t b = (truth & 2) ? ~0U : 0;
        const uint32_t c = (truth & 4) ? ~0U : 0;
        for (unsigned lane = 0; lane < VectorTraits<Value>::kLanes; ++lane)
        {
            input_a[lane] = a;
            input_b[lane] = b;
            input_c[lane] = c;
        }
        const Value result = Operation::Apply(
            VectorTraits<Value>::Load(input_a),
            VectorTraits<Value>::Load(input_b),
            VectorTraits<Value>::Load(input_c));
        VectorTraits<Value>::Store(output, result);
        for (unsigned lane = 0; lane < VectorTraits<Value>::kLanes; ++lane)
        {
            if (output[lane] != (a ^ b ^ c))
            {
                fprintf(stderr,
                    "Truth-table failure: width=%s operation=%s truth=%u lane=%u\n",
                    width, Operation::Name(), truth, lane);
                return false;
            }
        }
    }

    uint64_t random_state = 0x243f6a8885a308d3ULL;
    for (unsigned trial = 0; trial < 10000; ++trial)
    {
        for (unsigned lane = 0; lane < VectorTraits<Value>::kLanes; ++lane)
        {
            input_a[lane] = static_cast<uint32_t>(Mix64(random_state));
            input_b[lane] = static_cast<uint32_t>(Mix64(random_state));
            input_c[lane] = static_cast<uint32_t>(Mix64(random_state));
        }
        const Value result = Operation::Apply(
            VectorTraits<Value>::Load(input_a),
            VectorTraits<Value>::Load(input_b),
            VectorTraits<Value>::Load(input_c));
        VectorTraits<Value>::Store(output, result);
        for (unsigned lane = 0; lane < VectorTraits<Value>::kLanes; ++lane)
        {
            if (output[lane] != (input_a[lane] ^ input_b[lane] ^ input_c[lane]))
            {
                fprintf(stderr,
                    "Random failure: width=%s operation=%s trial=%u lane=%u\n",
                    width, Operation::Name(), trial, lane);
                return false;
            }
        }
    }
    return true;
}

template <typename Value>
static bool CheckWidth(const char* width)
{
    return CheckOperation<Value, ForcedXor2>(width) &&
        CheckOperation<Value, AutoXor2>(width) &&
        CheckOperation<Value, ExplicitXor3>(width) &&
        CheckOperation<Value, ProbeForcedXor2>(width) &&
        CheckOperation<Value, ProbeExplicitXor3>(width);
}

static uint8_t* AllocateAligned(uint64_t bytes)
{
#if defined(_WIN32)
    return static_cast<uint8_t*>(
        _aligned_malloc(static_cast<size_t>(bytes), 64));
#else
    void* pointer = 0;
    if (posix_memalign(&pointer, 64, static_cast<size_t>(bytes)) != 0)
        return 0;
    return static_cast<uint8_t*>(pointer);
#endif
}

static void FreeAligned(void* pointer)
{
#if defined(_WIN32)
    _aligned_free(pointer);
#else
    free(pointer);
#endif
}

static bool AllocateStream(StreamContext& context, uint64_t array_bytes)
{
    context.array_bytes = array_bytes;
    context.a = AllocateAligned(array_bytes);
    context.b = AllocateAligned(array_bytes);
    context.c = AllocateAligned(array_bytes);
    context.output = AllocateAligned(array_bytes);
    if (!context.a || !context.b || !context.c || !context.output)
    {
        FreeAligned(context.a);
        FreeAligned(context.b);
        FreeAligned(context.c);
        FreeAligned(context.output);
        return false;
    }

    uint64_t random_state = 0x13198a2e03707344ULL ^ array_bytes;
    for (uint64_t offset = 0; offset < array_bytes; ++offset)
    {
        context.a[offset] = static_cast<uint8_t>(Mix64(random_state));
        context.b[offset] = static_cast<uint8_t>(Mix64(random_state));
        context.c[offset] = static_cast<uint8_t>(Mix64(random_state));
        context.output[offset] = 0;
    }
    return true;
}

static void FreeStream(StreamContext& context)
{
    FreeAligned(context.a);
    FreeAligned(context.b);
    FreeAligned(context.c);
    FreeAligned(context.output);
}

template <typename Value>
static bool CheckStream(StreamContext& context)
{
    StreamKernel<Value, ForcedXor2>(1, &context);
    for (uint64_t offset = 0; offset < context.array_bytes; ++offset)
    {
        if (context.output[offset] !=
            static_cast<uint8_t>(context.a[offset] ^ context.b[offset] ^ context.c[offset]))
            return false;
    }
    StreamKernel<Value, AutoXor2>(1, &context);
    for (uint64_t offset = 0; offset < context.array_bytes; ++offset)
    {
        if (context.output[offset] !=
            static_cast<uint8_t>(context.a[offset] ^ context.b[offset] ^ context.c[offset]))
            return false;
    }
    StreamKernel<Value, ExplicitXor3>(1, &context);
    for (uint64_t offset = 0; offset < context.array_bytes; ++offset)
    {
        if (context.output[offset] !=
            static_cast<uint8_t>(context.a[offset] ^ context.b[offset] ^ context.c[offset]))
            return false;
    }
    return true;
}

template <typename Value>
static bool BenchmarkWidth(
    const char* width,
    bool quick,
    bool have_rdpru)
{
    const unsigned rounds = quick ? 3 : 9;
    const uint64_t latency_iterations = quick ? 8000000ULL : 32000000ULL;
    const uint64_t parallel_iterations = quick ? 2000000ULL : 8000000ULL;

#define LEO_MEASURE_DEPENDENCY(Shape, Iterations, Ops) \
    MeasurePair(width, #Shape, 0, \
        &Shape##Kernel<Value, ForcedXor2>, ForcedXor2::Name(), \
        &Shape##Kernel<Value, ExplicitXor3>, ExplicitXor3::Name(), \
        Iterations, Ops, 0, rounds, 0, have_rdpru); \
    MeasurePair(width, #Shape, 0, \
        &Shape##Kernel<Value, AutoXor2>, AutoXor2::Name(), \
        &Shape##Kernel<Value, ExplicitXor3>, ExplicitXor3::Name(), \
        Iterations, Ops, 0, rounds, 0, have_rdpru)

    LEO_MEASURE_DEPENDENCY(Latency, latency_iterations, 1.);
    LEO_MEASURE_DEPENDENCY(Parallel4, parallel_iterations, 4.);
#undef LEO_MEASURE_DEPENDENCY

    // Four arrays make the listed working-set size the aggregate footprint.
    static const uint64_t kWorkingSets[] = {
        4 * 1024ULL,
        32 * 1024ULL,
        256 * 1024ULL,
        1024 * 1024ULL,
        8 * 1024 * 1024ULL
    };
    const unsigned working_set_count = quick ? 3 :
        static_cast<unsigned>(sizeof(kWorkingSets) / sizeof(kWorkingSets[0]));
    const uint64_t target_logical_bytes = quick ?
        64 * 1024 * 1024ULL : 512 * 1024 * 1024ULL;
    for (unsigned i = 0; i < working_set_count; ++i)
    {
        StreamContext context;
        const uint64_t array_bytes = kWorkingSets[i] / 4;
        if (!AllocateStream(context, array_bytes))
        {
            fprintf(stderr, "Allocation failed for stream working set.\n");
            return false;
        }
        if (!CheckStream<Value>(context))
        {
            fprintf(stderr, "Streaming correctness failure for width %s.\n", width);
            FreeStream(context);
            return false;
        }

        uint64_t sweeps = target_logical_bytes / array_bytes;
        if (sweeps < 1)
            sweeps = 1;
        MeasurePair(width, "Stream", kWorkingSets[i],
            &StreamKernel<Value, ForcedXor2>, ForcedXor2::Name(),
            &StreamKernel<Value, ExplicitXor3>, ExplicitXor3::Name(),
            sweeps, static_cast<double>(array_bytes) / VectorTraits<Value>::kBytes,
            static_cast<double>(array_bytes), rounds, &context, have_rdpru);
        MeasurePair(width, "Stream", kWorkingSets[i],
            &StreamKernel<Value, AutoXor2>, AutoXor2::Name(),
            &StreamKernel<Value, ExplicitXor3>, ExplicitXor3::Name(),
            sweeps, static_cast<double>(array_bytes) / VectorTraits<Value>::kBytes,
            static_cast<double>(array_bytes), rounds, &context, have_rdpru);
        FreeStream(context);
    }
    return true;
}


} // namespace


int ff8xor_xor3_benchmark_run(
    bool quick,
    bool verify_only,
    unsigned pinned_cpu,
    bool have_rdpru)
{
    (void)pinned_cpu;
    if (!CheckWidth<__m128i>("XMM") ||
        !CheckWidth<__m256i>("YMM") ||
        !CheckWidth<__m512i>("ZMM"))
        return 1;

    StreamContext verification_context;
    if (!AllocateStream(verification_context, 4096))
        return 1;
    const bool stream_correct = CheckStream<__m128i>(verification_context) &&
        CheckStream<__m256i>(verification_context) &&
        CheckStream<__m512i>(verification_context);
    FreeStream(verification_context);
    if (!stream_correct)
    {
        fprintf(stderr, "Streaming correctness failure.\n");
        return 1;
    }

    printf("correctness=pass widths=%u operation_forms_per_width=%u "
           "truth_cases=%u random_vectors_per_form_width=%u "
           "stream_forms_widths=%u\n",
        3U, 5U, 8U, 10000U, 9U);
    if (verify_only)
        return 0;
    printf("NOTE: Stream ws is the aggregate four-array footprint; logical_GiB/s "
           "counts output bytes and modeled traffic is 4x.\n");
    printf("timing=steady_clock plus LFENCE/RDTSC reference samples\n");
    if (!BenchmarkWidth<__m128i>("XMM", quick, have_rdpru) ||
        !BenchmarkWidth<__m256i>("YMM", quick, have_rdpru) ||
        !BenchmarkWidth<__m512i>("ZMM", quick, have_rdpru))
        return 1;
    printf("sink=%llu\n", static_cast<unsigned long long>(BenchmarkSink));
    return 0;
}

#else


int ff8xor_xor3_benchmark_run(bool, bool, unsigned, bool)
{
    return 0;
}


#endif
