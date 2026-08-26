/*
    Selection-only direct-leaf screen for the register-resident T16 candidate.
    This is not promotion evidence: public dispatch, inactive neighbors, and
    same-binary end-to-end confidence gates are deliberately left to the later
    qualification campaign.
*/

#include "Leopard2Backend.h"
#include "LeopardFF8.h"
#include "leopard.h"

#include <algorithm>
#include <array>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <stdexcept>

namespace {

static const unsigned kSide = 16;
static const unsigned kShardBytes = 64;
static const size_t kPayloadBytes = kSide * kShardBytes;
static const unsigned kSamples = 31;
static const unsigned kWarmup = 4096;
static const unsigned kIterations = 32768;

void require(bool condition, const char* message)
{
    if (!condition)
        throw std::runtime_error(message);
}

uint64_t checksum(const uint8_t* data)
{
    uint64_t value = UINT64_C(0xcbf29ce484222325);
    for (size_t i = 0; i < kPayloadBytes; ++i)
        value = (value ^ data[i]) * UINT64_C(0x100000001b3);
    return value;
}

template<class Function>
double measure(Function function)
{
    const std::chrono::steady_clock::time_point begin =
        std::chrono::steady_clock::now();
    for (unsigned i = 0; i < kIterations; ++i)
        function();
    const std::chrono::steady_clock::time_point end =
        std::chrono::steady_clock::now();
    const double elapsed = std::chrono::duration_cast<
        std::chrono::duration<double, std::nano> >(end - begin).count();
    return elapsed / kIterations;
}

} // namespace

int main()
{
    try
    {
        require(leo_init() == Leopard_Success, "Leopard initialization failed");
        const leopard::backend::X86Features features =
            leopard::backend::DetectX86Features();
        require(features.avx512 && features.gfni,
            "AVX-512/GFNI is unavailable on this host");

        const leopard::backend::FF8K16R16B64PackedKernel candidate =
            leopard::backend::InitializeAVX512GFNIT16(
                leopard::ff8::MultiplyLogElement,
                leopard::ff8::SkewLogTable());
        require(candidate != NULL, "candidate KAT failed");

        alignas(64) std::array<uint8_t, kPayloadBytes> input;
        alignas(64) std::array<uint8_t, kPayloadBytes> baseline_output;
        alignas(64) std::array<uint8_t, kPayloadBytes> candidate_output;
        std::array<const void*, kSide> original;
        std::array<void*, kSide> recovery;
        uint64_t random = UINT64_C(0x7431362d7363726e);
        for (size_t i = 0; i < input.size(); ++i)
        {
            random += UINT64_C(0x9e3779b97f4a7c15);
            uint64_t value = random;
            value = (value ^ (value >> 30)) *
                UINT64_C(0xbf58476d1ce4e5b9);
            value = (value ^ (value >> 27)) *
                UINT64_C(0x94d049bb133111eb);
            value ^= value >> 31;
            input[i] = static_cast<uint8_t>(value);
        }
        for (unsigned row = 0; row < kSide; ++row)
        {
            original[row] = input.data() + row * kShardBytes;
            recovery[row] = baseline_output.data() + row * kShardBytes;
        }

        require(leopard::backend::TryAVX2FF8HighEncodeT16B64(
            original.data(), recovery.data()), "AVX2 baseline rejected input");
        candidate(input.data(), candidate_output.data());
        require(baseline_output == candidate_output,
            "candidate and AVX2 baseline outputs differ");

        for (unsigned i = 0; i < kWarmup; ++i)
        {
            require(leopard::backend::TryAVX2FF8HighEncodeT16B64(
                original.data(), recovery.data()), "AVX2 warmup failed");
            candidate(input.data(), candidate_output.data());
        }

        std::array<double, kSamples> baseline_ns;
        std::array<double, kSamples> candidate_ns;
        for (unsigned sample = 0; sample < kSamples; ++sample)
        {
            const auto baseline = [&]() {
                (void)leopard::backend::TryAVX2FF8HighEncodeT16B64(
                    original.data(), recovery.data());
            };
            const auto optimized = [&]() {
                candidate(input.data(), candidate_output.data());
            };
            if ((sample & 1U) == 0)
            {
                baseline_ns[sample] = measure(baseline);
                candidate_ns[sample] = measure(optimized);
            }
            else
            {
                candidate_ns[sample] = measure(optimized);
                baseline_ns[sample] = measure(baseline);
            }
        }

        std::array<double, kSamples> ratios;
        for (unsigned sample = 0; sample < kSamples; ++sample)
            ratios[sample] = baseline_ns[sample] / candidate_ns[sample];
        std::sort(baseline_ns.begin(), baseline_ns.end());
        std::sort(candidate_ns.begin(), candidate_ns.end());
        std::sort(ratios.begin(), ratios.end());
        const unsigned middle = kSamples / 2;
        const uint64_t baseline_digest = checksum(baseline_output.data());
        const uint64_t candidate_digest = checksum(candidate_output.data());
        std::printf(
            "{\"schema\":\"leopard2-t16-gfni-selection-v1\","
            "\"samples\":%u,\"iterations\":%u,"
            "\"baseline_median_ns\":%.9f,"
            "\"candidate_median_ns\":%.9f,"
            "\"median_paired_speedup\":%.9f,"
            "\"baseline_checksum\":%llu,"
            "\"candidate_checksum\":%llu}\n",
            kSamples, kIterations, baseline_ns[middle], candidate_ns[middle],
            ratios[middle],
            static_cast<unsigned long long>(baseline_digest),
            static_cast<unsigned long long>(candidate_digest));
        return 0;
    }
    catch (const std::exception& error)
    {
        std::fprintf(stderr, "%s\n", error.what());
        return 1;
    }
}
