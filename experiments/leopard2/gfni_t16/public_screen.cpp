/*
    Same-binary public-path selection screen for the exact T16 GFNI leaf.
    This remains diagnostic-only; the promotion runner separately freezes its
    binary, validates inactive neighbors, and computes paired confidence bounds.
*/

#include "Leopard2Direct.h"
#include "leopard2.h"

#include <algorithm>
#include <array>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <new>
#include <stdexcept>
#include <vector>

namespace {

static const unsigned kOriginalCount = 16;
static const unsigned kRecoveryCount = 16;
static const size_t kShardBytes = 64;
static const size_t kPayloadBytes = kOriginalCount * kShardBytes;
static const unsigned kSamples = 31;
static const unsigned kWarmup = 256;
static const unsigned kIterations = 8192;

void require(bool condition, const char* message)
{
    if (!condition)
        throw std::runtime_error(message);
}

void require_result(leo2_result actual, leo2_result expected, const char* message)
{
    if (actual != expected)
        throw std::runtime_error(message);
}

class AlignedBuffer
{
public:
    explicit AlignedBuffer(size_t bytes)
        : pointer_(NULL)
        , bytes_(bytes)
    {
        if (posix_memalign(&pointer_, leo2_scratch_alignment(), bytes) != 0)
            pointer_ = NULL;
        if (!pointer_)
            throw std::bad_alloc();
        std::memset(pointer_, 0, bytes);
    }

    ~AlignedBuffer() { std::free(pointer_); }
    void* data() const { return pointer_; }
    size_t size() const { return bytes_; }

private:
    AlignedBuffer(const AlignedBuffer&);
    AlignedBuffer& operator=(const AlignedBuffer&);
    void* pointer_;
    size_t bytes_;
};

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
    return std::chrono::duration_cast<
        std::chrono::duration<double, std::nano> >(end - begin).count() /
        kIterations;
}

void normalize_mode(bool enabled)
{
    require(leopard2_internal::
            SetK16R16B64AVX512GFNIEnabledForDiagnostics(enabled),
        "failed to set T16 GFNI diagnostic mode");
    require(leopard2_internal::
            FinishK16R16B64AVX512GFNIRouteProbeForDiagnostics(),
        "failed to normalize T16 GFNI diagnostic mode");
}

struct Summary
{
    double control_ns;
    double candidate_ns;
    double speedup;
};

template<class Function>
Summary screen(Function function)
{
    for (unsigned i = 0; i < kWarmup; ++i)
    {
        normalize_mode(false);
        function();
        normalize_mode(true);
        function();
    }

    std::array<double, kSamples> control;
    std::array<double, kSamples> candidate;
    std::array<double, kSamples> ratio;
    for (unsigned sample = 0; sample < kSamples; ++sample)
    {
        if ((sample & 1U) == 0)
        {
            normalize_mode(false);
            control[sample] = measure(function);
            normalize_mode(true);
            candidate[sample] = measure(function);
        }
        else
        {
            normalize_mode(true);
            candidate[sample] = measure(function);
            normalize_mode(false);
            control[sample] = measure(function);
        }
        ratio[sample] = control[sample] / candidate[sample];
    }
    std::sort(control.begin(), control.end());
    std::sort(candidate.begin(), candidate.end());
    std::sort(ratio.begin(), ratio.end());
    const unsigned middle = kSamples / 2;
    Summary result = { control[middle], candidate[middle], ratio[middle] };
    return result;
}

} // namespace

int main()
{
    try
    {
        leo2_context_options options = {};
        options.struct_size = sizeof(options);
        options.backend = LEO2_BACKEND_AUTO;
        options.thread_count = 1;
        leo2_context* context = NULL;
        require_result(leo2_context_create(&options, &context), LEO2_SUCCESS,
            "failed to create AUTO context");
        leo2_codec* codec = NULL;
        require_result(leo2_codec_create(context, kOriginalCount,
            kRecoveryCount, LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8,
            NULL, &codec), LEO2_SUCCESS, "failed to create T16 codec");
        require(leopard2_internal::
                K16R16B64AVX512GFNIAvailableForDiagnostics(codec),
            "T16 GFNI candidate is unavailable");

        size_t scratch_bytes = 0;
        require_result(leo2_encode_scratch_size(codec, kShardBytes,
            &scratch_bytes), LEO2_SUCCESS, "failed to query scratch");
        AlignedBuffer scratch(scratch_bytes);
        std::vector<uint8_t> input_storage(kPayloadBytes + 8, 0);
        std::vector<uint8_t> output_storage(kPayloadBytes + 8, 0);
        uint8_t* const input = &input_storage[1];
        uint8_t* const output = &output_storage[3];
        std::array<const void*, kOriginalCount> original;
        std::array<void*, kRecoveryCount> recovery;
        uint64_t random = UINT64_C(0x7431362d7075626c);
        for (size_t i = 0; i < kPayloadBytes; ++i)
        {
            random ^= random << 13;
            random ^= random >> 7;
            random ^= random << 17;
            input[i] = static_cast<uint8_t>(random >> 24);
        }
        for (unsigned row = 0; row < kOriginalCount; ++row)
        {
            original[row] = input + row * kShardBytes;
            recovery[row] = output + row * kShardBytes;
        }

        const auto ordinary = [&]() {
            const leo2_result result = leo2_encode(codec, kShardBytes,
                original.data(), recovery.data(), scratch.data(), scratch.size());
            if (result != LEO2_SUCCESS)
                std::abort();
        };
        leo2_encode_batch_item item = {
            kShardBytes, original.data(), recovery.data(),
            scratch.data(), scratch.size()
        };
        const auto batch_one = [&]() {
            if (leo2_encode_batch(codec, &item, 1) != LEO2_SUCCESS)
                std::abort();
        };

        normalize_mode(false);
        ordinary();
        const uint64_t control_digest = checksum(output);
        normalize_mode(true);
        ordinary();
        const uint64_t candidate_digest = checksum(output);
        require(control_digest == candidate_digest,
            "public control and candidate checksums differ");
        require(leopard2_internal::K16R16B64AVX512GFNISelectedForDiagnostics(
                codec, kShardBytes), "candidate selector is not active");

        const Summary ordinary_summary = screen(ordinary);
        const Summary batch_summary = screen(batch_one);
        normalize_mode(false);
        std::printf(
            "{\"schema\":\"leopard2-t16-gfni-public-selection-v1\","
            "\"samples\":%u,\"iterations\":%u,"
            "\"ordinary\":{\"control_ns\":%.9f,"
            "\"candidate_ns\":%.9f,\"speedup\":%.9f},"
            "\"batch_one\":{\"control_ns\":%.9f,"
            "\"candidate_ns\":%.9f,\"speedup\":%.9f},"
            "\"control_checksum\":%llu,\"candidate_checksum\":%llu}\n",
            kSamples, kIterations,
            ordinary_summary.control_ns, ordinary_summary.candidate_ns,
            ordinary_summary.speedup,
            batch_summary.control_ns, batch_summary.candidate_ns,
            batch_summary.speedup,
            static_cast<unsigned long long>(control_digest),
            static_cast<unsigned long long>(candidate_digest));
        leo2_codec_destroy(codec);
        leo2_context_destroy(context);
        return 0;
    }
    catch (const std::exception& error)
    {
        std::fprintf(stderr, "%s\n", error.what());
        return 1;
    }
}
