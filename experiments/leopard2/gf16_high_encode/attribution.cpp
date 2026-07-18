/*
    Diagnostic attribution harness for the legacy-high GF16 encoder.

    This deliberately compares the public Leopard2 validation/pointer setup
    with the exact same production transform called through its private field
    entry point.  It is not a throughput benchmark or public API example.
*/

#include "Leopard2Backend.h"
#include "LeopardFF16.h"
#include "leopard.h"
#include "leopard2.h"

#include <algorithm>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <stdint.h>
#include <vector>

namespace {

class Buffer
{
public:
    Buffer() : data_(NULL), bytes_(0) {}
    ~Buffer() { free(data_); }

    bool Reset(size_t bytes)
    {
        free(data_);
        data_ = NULL;
        bytes_ = bytes;
        if (posix_memalign(&data_, 64, bytes) != 0)
            return false;
        memset(data_, 0, bytes);
        return true;
    }

    uint8_t* bytes() { return static_cast<uint8_t*>(data_); }
    void* data() { return data_; }
    size_t size() const { return bytes_; }

private:
    Buffer(const Buffer&);
    Buffer& operator=(const Buffer&);
    void* data_;
    size_t bytes_;
};

uint32_t NextPow2(uint32_t value)
{
    uint32_t result = 1;
    while (result < value)
        result <<= 1;
    return result;
}

double Median(std::vector<double> values)
{
    std::sort(values.begin(), values.end());
    return values[values.size() / 2];
}

template<class Callable>
double Measure(unsigned iterations, unsigned reuse, Callable callable)
{
    typedef std::chrono::steady_clock Clock;
    std::vector<double> samples;
    samples.reserve(iterations);
    for (unsigned sample = 0; sample < iterations; ++sample)
    {
        const Clock::time_point begin = Clock::now();
        for (unsigned i = 0; i < reuse; ++i)
            callable();
        const Clock::time_point end = Clock::now();
        samples.push_back(std::chrono::duration_cast<
            std::chrono::duration<double, std::micro> >(end - begin).count()
            / reuse);
    }
    return Median(samples);
}

void Fail(const char* message)
{
    fprintf(stderr, "%s\n", message);
    exit(1);
}

leo2_backend ParseBackend(const char* text)
{
    if (strcmp(text, "auto") == 0)
        return LEO2_BACKEND_AUTO;
    if (strcmp(text, "avx2") == 0)
        return LEO2_BACKEND_AVX2;
    if (strcmp(text, "avx512") == 0 || strcmp(text, "avx512vl") == 0)
        return LEO2_BACKEND_AVX512;
    Fail("backend must be auto, avx2, or avx512");
    return LEO2_BACKEND_AUTO;
}

} // namespace

int main(int argc, char** argv)
{
    if (argc != 4 && argc != 5)
    {
        fprintf(stderr,
            "usage: %s K R shard_bytes [auto|avx2|avx512]\n", argv[0]);
        return 2;
    }
    const uint32_t k = static_cast<uint32_t>(strtoul(argv[1], NULL, 10));
    const uint32_t r = static_cast<uint32_t>(strtoul(argv[2], NULL, 10));
    const size_t shard_bytes = static_cast<size_t>(strtoull(argv[3], NULL, 10));
    const uint32_t side = NextPow2(r);
    const leo2_backend requested_backend = argc == 5
        ? ParseBackend(argv[4]) : LEO2_BACKEND_AUTO;

    leo2_context_options context_options = {};
    context_options.struct_size = sizeof(context_options);
    context_options.backend = requested_backend;
    context_options.thread_count = 1;
    leo2_context* context = NULL;
    if (leo2_context_create(&context_options, &context) != LEO2_SUCCESS)
        Fail("context creation failed");
    leo2_codec_options codec_options = {};
    codec_options.struct_size = sizeof(codec_options);
    leo2_codec* codec = NULL;
    if (leo2_codec_create(context, k, r, LEO2_PROFILE_LEGACY_HIGH_V1,
            LEO2_FIELD_GF16, &codec_options, &codec) != LEO2_SUCCESS)
        Fail("codec creation failed");

    size_t scratch_bytes = 0;
    if (leo2_encode_scratch_size(codec, shard_bytes, &scratch_bytes) !=
            LEO2_SUCCESS)
        Fail("scratch query failed");

    Buffer originals;
    Buffer public_outputs;
    Buffer public_scratch;
    Buffer direct_outputs;
    Buffer direct_work;
    Buffer contiguous_work;
    if (!originals.Reset(static_cast<size_t>(k) * shard_bytes) ||
        !public_outputs.Reset(static_cast<size_t>(r) * shard_bytes) ||
        !public_scratch.Reset(scratch_bytes) ||
        !direct_outputs.Reset(static_cast<size_t>(r) * shard_bytes) ||
        !direct_work.Reset(static_cast<size_t>(side) * 2 * shard_bytes) ||
        !contiguous_work.Reset(static_cast<size_t>(side) * 2 * shard_bytes))
        Fail("allocation failed");

    uint64_t random = UINT64_C(0x9e3779b97f4a7c15);
    for (size_t i = 0; i < originals.size(); ++i)
    {
        random ^= random << 13;
        random ^= random >> 7;
        random ^= random << 17;
        originals.bytes()[i] = static_cast<uint8_t>(random >> 56);
    }

    std::vector<const void*> original(k);
    std::vector<void*> public_recovery(r);
    std::vector<void*> direct_split(side * 2);
    std::vector<void*> direct_contiguous(side * 2);
    for (uint32_t i = 0; i < k; ++i)
        original[i] = originals.bytes() + static_cast<size_t>(i) * shard_bytes;
    for (uint32_t i = 0; i < r; ++i)
    {
        public_recovery[i] = public_outputs.bytes() +
            static_cast<size_t>(i) * shard_bytes;
        direct_split[i] = direct_outputs.bytes() +
            static_cast<size_t>(i) * shard_bytes;
    }
    for (uint32_t i = r; i < side * 2; ++i)
        direct_split[i] = direct_work.bytes() +
            static_cast<size_t>(i) * shard_bytes;
    for (uint32_t i = 0; i < side * 2; ++i)
        direct_contiguous[i] = contiguous_work.bytes() +
            static_cast<size_t>(i) * shard_bytes;

    leo2_encode_batch_item item = {};
    item.shard_bytes = shard_bytes;
    item.original = &original[0];
    item.recovery = &public_recovery[0];
    item.scratch = public_scratch.data();
    item.scratch_bytes = public_scratch.size();

    const leopard::backend::Ops* const selected_ops = requested_backend ==
            LEO2_BACKEND_AUTO
        ? &leopard::backend::GetDefaultOps()
        : leopard::backend::GetQualifiedOps(requested_backend);
    if (!selected_ops)
        Fail("selected backend ops are unavailable");
    const leopard::backend::Ops& ops = *selected_ops;
    const unsigned warmup = 2;
    for (unsigned i = 0; i < warmup; ++i)
    {
        if (leo2_encode(codec, shard_bytes, &original[0],
                &public_recovery[0], public_scratch.data(),
                public_scratch.size()) != LEO2_SUCCESS)
            Fail("public encode failed");
        leopard::ff16::ReedSolomonEncode(ops, shard_bytes, k, r, r, side,
            &original[0], &direct_split[0], NULL);
        leopard::ff16::ReedSolomonEncode(ops, shard_bytes, k, r, r, side,
            &original[0], &direct_contiguous[0], NULL);
    }
    if (memcmp(public_outputs.data(), direct_outputs.data(),
            static_cast<size_t>(r) * shard_bytes) != 0 ||
        memcmp(public_outputs.data(), contiguous_work.data(),
            static_cast<size_t>(r) * shard_bytes) != 0)
        Fail("parity mismatch");

    const unsigned iterations = 9;
    const unsigned reuse = 2;
    const double public_one = Measure(iterations, reuse, [&]() {
        if (leo2_encode(codec, shard_bytes, &original[0],
                &public_recovery[0], public_scratch.data(),
                public_scratch.size()) != LEO2_SUCCESS)
            Fail("public encode failed");
    });
    const double public_batch = Measure(iterations, reuse, [&]() {
        if (leo2_encode_batch(codec, &item, 1) != LEO2_SUCCESS)
            Fail("public batch encode failed");
    });
    const double split = Measure(iterations, reuse, [&]() {
        leopard::ff16::ReedSolomonEncode(ops, shard_bytes, k, r, r, side,
            &original[0], &direct_split[0], NULL);
    });
    const double contiguous = Measure(iterations, reuse, [&]() {
        leopard::ff16::ReedSolomonEncode(ops, shard_bytes, k, r, r, side,
            &original[0], &direct_contiguous[0], NULL);
    });

    printf("K=%u R=%u bytes=%zu backend=%s public_one_us=%.3f "
           "public_batch_us=%.3f split_us=%.3f contiguous_us=%.3f "
           "split_over_public=%.6f contiguous_over_split=%.6f\n",
        k, r, shard_bytes, ops.name, public_one, public_batch, split,
        contiguous, split / public_one, contiguous / split);

    leo2_codec_destroy(codec);
    leo2_context_destroy(context);
    return 0;
}
