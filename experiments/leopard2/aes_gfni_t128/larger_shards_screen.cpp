/*
    Selection-only screen for the larger-shard K=65/R=65/T=128 native-Cantor
    AVX-512/GFNI candidate.  These timings are not promotion evidence; the
    retained benchmark schema and fresh confirmation must be built after the
    implementation policy is frozen.
*/

#include "Leopard2Direct.h"
#include "leopard2.h"

#include <algorithm>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <memory>
#include <new>
#include <stdexcept>
#include <vector>

namespace {

static const unsigned kOriginalCount = 65;
static const unsigned kRecoveryCount = 65;
static volatile uint64_t Sink = 0;

class AlignedBuffer
{
public:
    explicit AlignedBuffer(size_t bytes)
        : data_(NULL), bytes_(bytes)
    {
        if (posix_memalign(&data_, 64, bytes_) != 0 || !data_)
            throw std::bad_alloc();
    }

    ~AlignedBuffer()
    {
        std::free(data_);
    }

    uint8_t* bytes()
    {
        return static_cast<uint8_t*>(data_);
    }

    size_t size() const
    {
        return bytes_;
    }

private:
    void* data_;
    size_t bytes_;
};

void Require(bool condition, const char* message)
{
    if (!condition)
        throw std::runtime_error(message);
}

void RequireSuccess(leo2_result result, const char* message)
{
    if (result != LEO2_SUCCESS)
        throw std::runtime_error(message);
}

void Fill(uint8_t* data, size_t bytes, uint64_t seed)
{
    uint64_t state = seed;
    for (size_t i = 0; i < bytes; ++i)
    {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        data[i] = static_cast<uint8_t>(state >> 24);
    }
}

double Median(std::vector<double> values)
{
    std::sort(values.begin(), values.end());
    const size_t middle = values.size() / 2;
    return values.size() & 1U
        ? values[middle]
        : (values[middle - 1] + values[middle]) * 0.5;
}

class Workload
{
public:
    explicit Workload(size_t shard_bytes)
        : shard_bytes_(shard_bytes)
        , input_(kOriginalCount * shard_bytes)
        , output_(kRecoveryCount * shard_bytes)
        , context_(NULL)
        , codec_(NULL)
    {
        leo2_context_options options = {};
        options.struct_size = sizeof(options);
        options.backend = LEO2_BACKEND_AUTO;
        options.thread_count = 1;
        RequireSuccess(leo2_context_create(&options, &context_),
            "create AUTO context");
        RequireSuccess(leo2_codec_create(context_,
            kOriginalCount, kRecoveryCount,
            LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8,
            NULL, &codec_), "create K65/R65 codec");
        Require(leopard2_internal::
                SetK65R65T128AVX512GFNIEnabledForDiagnostics(true),
            "enable larger GFNI selector");
        Require(leopard2_internal::
                K65R65T128AVX512GFNILargerSelectedForDiagnostics(
                    codec_, shard_bytes_),
            "larger GFNI selector is unavailable");
        Require(leopard2_internal::
                FinishK65R65T128AVX512GFNIRouteProbeForDiagnostics(),
            "finish selector availability probe");

        size_t scratch_bytes = 0;
        RequireSuccess(leo2_encode_scratch_size(
            codec_, shard_bytes_, &scratch_bytes), "query scratch");
        scratch_.reset(new AlignedBuffer(scratch_bytes));
        control_.resize(kRecoveryCount * shard_bytes_);
        original_.resize(kOriginalCount);
        recovery_.resize(kRecoveryCount);
        for (unsigned i = 0; i < kOriginalCount; ++i)
        {
            original_[i] = input_.bytes() +
                static_cast<size_t>(i) * shard_bytes_;
        }
        for (unsigned i = 0; i < kRecoveryCount; ++i)
        {
            recovery_[i] = output_.bytes() +
                static_cast<size_t>(i) * shard_bytes_;
        }
        item_.shard_bytes = shard_bytes_;
        item_.original = &original_[0];
        item_.recovery = &recovery_[0];
        item_.scratch = scratch_->bytes();
        item_.scratch_bytes = scratch_->size();
        Fill(input_.bytes(), kOriginalCount * shard_bytes_,
            UINT64_C(0x6c61726765723634) + shard_bytes_);
        CheckCorrectness();
    }

    ~Workload()
    {
        leo2_codec_destroy(codec_);
        leo2_context_destroy(context_);
    }

    double Measure(bool candidate, bool batch, unsigned reuse)
    {
        Require(leopard2_internal::
                SetK65R65T128AVX512GFNIEnabledForDiagnostics(candidate),
            "set diagnostic mode");
        Execute(batch);
        const unsigned calls = leopard2_internal::
            K65R65T128AVX512GFNICallCountForDiagnostics();
        const size_t tiles = leopard2_internal::
            K65R65T128AVX512GFNITileCountForDiagnostics();
        Require(calls == (candidate ? 1U : 0U),
            "route-probe call count mismatch");
        Require(tiles == (candidate ? shard_bytes_ / 64U : 0U),
            "route-probe tile count mismatch");
        Require(leopard2_internal::
                FinishK65R65T128AVX512GFNIRouteProbeForDiagnostics(),
            "finish route probe");

        for (unsigned i = 0; i < 32; ++i)
            Execute(batch);
        const std::chrono::steady_clock::time_point started =
            std::chrono::steady_clock::now();
        for (unsigned i = 0; i < reuse; ++i)
            Execute(batch);
        const std::chrono::steady_clock::time_point stopped =
            std::chrono::steady_clock::now();
        Sink ^= output_.bytes()[
            (static_cast<size_t>(reuse) * 131U) % output_.size()];
        return std::chrono::duration<double, std::nano>(
            stopped - started).count() / reuse;
    }

private:
    void Execute(bool batch)
    {
        const leo2_result result = batch
            ? leo2_encode_batch(codec_, &item_, 1)
            : leo2_encode(codec_, shard_bytes_, &original_[0], &recovery_[0],
                scratch_->bytes(), scratch_->size());
        RequireSuccess(result, "execute encode");
    }

    void CheckCorrectness()
    {
        Require(leopard2_internal::
                SetK65R65T128AVX512GFNIEnabledForDiagnostics(false),
            "disable candidate for correctness");
        Execute(false);
        std::memcpy(&control_[0], output_.bytes(), control_.size());
        Require(leopard2_internal::
                FinishK65R65T128AVX512GFNIRouteProbeForDiagnostics(),
            "finish control correctness probe");
        Require(leopard2_internal::
                SetK65R65T128AVX512GFNIEnabledForDiagnostics(true),
            "enable candidate for correctness");
        Execute(false);
        Require(std::memcmp(&control_[0], output_.bytes(),
                control_.size()) == 0,
            "candidate output differs from mature control");
        Require(leopard2_internal::
                FinishK65R65T128AVX512GFNIRouteProbeForDiagnostics(),
            "finish candidate correctness probe");
    }

    size_t shard_bytes_;
    AlignedBuffer input_;
    AlignedBuffer output_;
    std::unique_ptr<AlignedBuffer> scratch_;
    std::vector<uint8_t> control_;
    std::vector<const void*> original_;
    std::vector<void*> recovery_;
    leo2_encode_batch_item item_;
    leo2_context* context_;
    leo2_codec* codec_;
};

void Screen(size_t shard_bytes, bool batch)
{
    Workload workload(shard_bytes);
    const unsigned reuse = static_cast<unsigned>(524288U / shard_bytes);
    std::vector<double> control;
    std::vector<double> candidate;
    for (unsigned round = 0; round < 9; ++round)
    {
        const bool reverse = (round & 1U) != 0;
        for (unsigned index = 0; index < 4; ++index)
        {
            const bool use_candidate = reverse
                ? (index == 0 || index == 3)
                : (index == 1 || index == 2);
            const double elapsed = workload.Measure(
                use_candidate, batch, reuse);
            (use_candidate ? candidate : control).push_back(elapsed);
        }
    }
    const double control_ns = Median(control);
    const double candidate_ns = Median(candidate);
    std::printf(
        "B=%zu api=%s reuse=%u control_ns=%.6f candidate_ns=%.6f "
        "speedup=%.6fx\n",
        shard_bytes, batch ? "batch1" : "one_shot", reuse,
        control_ns, candidate_ns, control_ns / candidate_ns);
}

} // namespace

int main()
{
    try
    {
        static const size_t kBytes[] = {
            128, 256, 512, 1024, 2048, 4096
        };
        for (size_t i = 0; i < sizeof(kBytes) / sizeof(kBytes[0]); ++i)
        {
            Screen(kBytes[i], true);
            Screen(kBytes[i], false);
        }
        std::printf("sink=%llu\n",
            static_cast<unsigned long long>(Sink));
        return 0;
    }
    catch (const std::exception& error)
    {
        std::fprintf(stderr, "larger-shard screen failure: %s\n", error.what());
        return 1;
    }
}
