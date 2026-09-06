// Bounded same-source diagnostic for leopard-79h.38.5.4.9.
// Not an exact-Leopard1 benchmark or a production qualification harness.
#include "leopard2.h"
#include "Leopard2Direct.h"

#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <memory>
#include <stdexcept>
#include <vector>

namespace {
void Require(bool value, const char* what)
{
    if (!value) throw std::runtime_error(what);
}
void Check(leo2_result result) { Require(result == LEO2_SUCCESS, "API failure"); }

struct Aligned
{
    void* data;
    explicit Aligned(size_t bytes) : data(NULL)
    {
        Require(posix_memalign(&data, 64, bytes ? bytes : 64) == 0, "allocation");
        std::memset(data, 0, bytes);
    }
    ~Aligned() { std::free(data); }
    Aligned(const Aligned&) = delete;
    Aligned& operator=(const Aligned&) = delete;
};

uint64_t Hash(const std::vector<uint8_t>& data)
{
    uint64_t hash = UINT64_C(14695981039346656037);
    for (uint8_t value : data)
        hash = (hash ^ value) * UINT64_C(1099511628211);
    return hash;
}

struct Cell { unsigned k, r; size_t bytes; leo2_backend backend; };
const Cell kCells[] = {
    {1000, 200, 32768, LEO2_BACKEND_AVX512},
    {1000, 200, 65536, LEO2_BACKEND_AVX512},
    {1000, 199, 65536, LEO2_BACKEND_AVX512},
    {4096, 512, 4096, LEO2_BACKEND_AVX512},
    {1000, 200, 65536, LEO2_BACKEND_AVX2},
    {1000, 200, 65536, LEO2_BACKEND_AUTO}
};
}

int main(int argc, char** argv)
{
    try
    {
        Require(argc == 3 && (std::strcmp(argv[1], "--check") == 0 ||
            std::strcmp(argv[1], "--measure") == 0),
            "usage: split_cache_screen --check|--measure cell[0..5]");
        Require(std::strlen(argv[2]) == 1 && argv[2][0] >= '0' &&
            argv[2][0] <= '5', "invalid cell");
        const unsigned index = static_cast<unsigned>(argv[2][0] - '0');
        const Cell& cell = kCells[index];
        const bool measure = std::strcmp(argv[1], "--measure") == 0;
        leo2_context_options options = {};
        options.struct_size = sizeof(options);
        options.backend = cell.backend;
        options.thread_count = 1;
        leo2_context* raw_context = NULL;
        Check(leo2_context_create(&options, &raw_context));
        std::unique_ptr<leo2_context, decltype(&leo2_context_destroy)>
            context(raw_context, leo2_context_destroy);
        leo2_codec* raw_codec = NULL;
        Check(leo2_codec_create(context.get(), cell.k, cell.r,
            LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF16, NULL, &raw_codec));
        std::unique_ptr<leo2_codec, decltype(&leo2_codec_destroy)>
            codec(raw_codec, leo2_codec_destroy);
        const bool gfni = leopard2_internal::
            AutoGF16GFNIEncodeSelectedForDiagnostics(codec.get(), cell.bytes);
        Require(gfni == (index == 5), "unexpected AUTO GFNI selection");
        const leo2_backend reported = leo2_context_backend(context.get());
        Require(reported == (index < 4 ? LEO2_BACKEND_AVX512 :
            LEO2_BACKEND_AVX2), "unexpected context backend");
        const char* route = gfni ? "gfni" : index < 4 ? "avx512" : "avx2";
        std::vector<uint8_t> source(static_cast<size_t>(cell.k) * cell.bytes);
        std::vector<uint8_t> output(static_cast<size_t>(cell.r) * cell.bytes);
        std::vector<const void*> inputs(cell.k);
        std::vector<void*> outputs(cell.r);
        uint32_t random = 20260906;
        for (uint8_t& value : source)
        {
            random ^= random << 13;
            random ^= random >> 17;
            random ^= random << 5;
            value = static_cast<uint8_t>(random);
        }
        for (unsigned i = 0; i < cell.k; ++i)
            inputs[i] = source.data() + static_cast<size_t>(i) * cell.bytes;
        for (unsigned i = 0; i < cell.r; ++i)
            outputs[i] = output.data() + static_cast<size_t>(i) * cell.bytes;
        size_t scratch_bytes = 0;
        Check(leo2_encode_scratch_size(codec.get(), cell.bytes, &scratch_bytes));
        Aligned scratch(scratch_bytes);
        const uint64_t input_hash = Hash(source);
        const auto encode = [&]() {
            Check(leo2_encode(codec.get(), cell.bytes, inputs.data(),
                outputs.data(), scratch.data, scratch_bytes));
        };
        encode();
        const uint64_t output_hash = Hash(output);
        // Timing excludes construction, allocation and digest calculation.
        // Check mode never reads a clock and emits no performance samples.
        std::vector<long long> samples;
        if (measure)
        {
            for (unsigned i = 0; i < 4; ++i) encode();
            samples.reserve(21);
            for (unsigned i = 0; i < 21; ++i)
            {
                const auto start = std::chrono::steady_clock::now();
                encode();
                const auto end = std::chrono::steady_clock::now();
                samples.push_back(std::chrono::duration_cast<
                    std::chrono::nanoseconds>(end - start).count());
            }
        }
        Require(Hash(source) == input_hash && Hash(output) == output_hash,
            "workload changed");
        Require(leopard2_internal::AutoGF16GFNIEncodeSelectedForDiagnostics(
            codec.get(), cell.bytes) == gfni, "route changed");
        std::printf("{\"schema\":\"leopard2-gf16-split-screen/v1\","
            "\"cell\":%u,\"k\":%u,\"r\":%u,\"bytes\":%zu,"
            "\"execution_route\":\"%s\",\"scratch_bytes\":%zu,"
            "\"input_hash\":\"%016llx\",\"output_hash\":\"%016llx\","
            "\"samples_ns\":[", index, cell.k, cell.r, cell.bytes,
            route, scratch_bytes, static_cast<unsigned long long>(input_hash),
            static_cast<unsigned long long>(output_hash));
        for (size_t i = 0; i < samples.size(); ++i)
            std::printf("%s%lld", i ? "," : "", samples[i]);
        std::puts("]}");
        return 0;
    }
    catch (const std::exception& error)
    {
        std::fprintf(stderr, "%s\n", error.what());
        return 1;
    }
}
