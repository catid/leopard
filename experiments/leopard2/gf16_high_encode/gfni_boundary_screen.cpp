// leopard-79h.38.5.4.17: unchanged-codec, explicit-backend boundary diagnostic.
#include "leopard2.h"
#include "Leopard2Direct.h"
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <memory>
#include <stdexcept>
#include <vector>

#ifndef LEO_BOUNDARY_CODEC_COMMIT
#error Supply the independently pinned codec identity.
#endif

namespace {
void Require(bool ok, const char* message)
{
    if (!ok) throw std::runtime_error(message);
}

struct Aligned
{
    uint8_t* data;
    explicit Aligned(size_t bytes) : data(NULL)
    {
        void* allocation = NULL;
        Require(posix_memalign(&allocation, 64, bytes ? bytes : 64) == 0, "allocation");
        data = static_cast<uint8_t*>(allocation);
        std::memset(data, 0, bytes);
    }
    ~Aligned() { std::free(data); }
    Aligned(const Aligned&) = delete;
    Aligned& operator=(const Aligned&) = delete;
};

uint64_t Hash(const uint8_t* data, size_t size)
{
    uint64_t result = UINT64_C(14695981039346656037);
    for (size_t i = 0; i < size; ++i)
        result = (result ^ data[i]) * UINT64_C(1099511628211);
    return result;
}

void Fill(uint8_t* data, size_t size, uint32_t& random)
{
    for (size_t i = 0; i < size; ++i)
    {
        random ^= random << 13;
        random ^= random >> 17;
        random ^= random << 5;
        data[i] = static_cast<uint8_t>(random);
    }
}

struct Codec
{
    typedef std::unique_ptr<leo2_context, decltype(&leo2_context_destroy)> Context;
    typedef std::unique_ptr<leo2_codec, decltype(&leo2_codec_destroy)> Handle;
    Context context;
    Handle codec;
    Codec(unsigned k, unsigned r, leo2_field field, leo2_backend backend)
        : context(NULL, leo2_context_destroy), codec(NULL, leo2_codec_destroy)
    {
        leo2_context_options options = {};
        options.struct_size = sizeof(options);
        options.backend = backend;
        options.thread_count = 1;
        leo2_context* raw_context = NULL;
        Require(leo2_context_create(&options, &raw_context) == LEO2_SUCCESS, "context");
        context.reset(raw_context);
        Require(leo2_context_backend(context.get()) ==
            (backend == LEO2_BACKEND_AUTO ? LEO2_BACKEND_AVX2 : backend), "backend");
        leo2_codec* raw_codec = NULL;
        Require(leo2_codec_create(context.get(), k, r, LEO2_PROFILE_LEGACY_HIGH_V1,
            field, NULL, &raw_codec) == LEO2_SUCCESS, "codec");
        codec.reset(raw_codec);
    }
    size_t Scratch(size_t bytes) const
    {
        size_t result = 0;
        Require(leo2_encode_scratch_size(codec.get(), bytes, &result) == LEO2_SUCCESS,
            "scratch query");
        return result;
    }
};

struct Cell { unsigned k, r; size_t bytes; };
const Cell kCells[] = {{1000, 200, 32768}, {1000, 199, 65536}};
}

#ifndef LEO_BOUNDARY_NO_MAIN
int main(int argc, char** argv)
{
    try
    {
        Require((argc == 4 || argc == 5) &&
            (!std::strcmp(argv[1], "--check") || !std::strcmp(argv[1], "--measure")),
            "usage: --check|--measure cell[0..1] avx2|gfni|auto [untimed_parity_file]");
        Require(std::strlen(argv[2]) == 1 && argv[2][0] >= '0' && argv[2][0] <= '1',
            "cell");
        const unsigned index = static_cast<unsigned>(argv[2][0] - '0');
        const Cell cell = kCells[index];
        const bool measured = !std::strcmp(argv[1], "--measure");
        Require(!measured || argc == 4, "no timed parity dump");
        Require(!std::strcmp(argv[3], "avx2") || !std::strcmp(argv[3], "gfni") ||
            !std::strcmp(argv[3], "auto"), "requested backend");
        const leo2_backend backend = !std::strcmp(argv[3], "gfni") ? LEO2_BACKEND_GFNI :
            !std::strcmp(argv[3], "auto") ? LEO2_BACKEND_AUTO : LEO2_BACKEND_AVX2;
        Codec codec(cell.k, cell.r, LEO2_FIELD_GF16, backend);
        Require(!leopard2_internal::AutoGF16GFNIEncodeSelectedForDiagnostics(
            codec.codec.get(), cell.bytes), "AUTO boundary unexpectedly widened");
        const size_t input_bytes = static_cast<size_t>(cell.k) * cell.bytes;
        const size_t output_bytes = static_cast<size_t>(cell.r) * cell.bytes;
        const size_t scratch_bytes = codec.Scratch(cell.bytes);
        Aligned source(input_bytes), output(output_bytes), scratch(scratch_bytes);
        uint32_t random = 20260906;
        Fill(source.data, input_bytes, random);
        std::vector<const void*> inputs(cell.k);
        std::vector<void*> outputs(cell.r);
        for (unsigned i = 0; i < cell.k; ++i) inputs[i] = source.data + i * cell.bytes;
        for (unsigned i = 0; i < cell.r; ++i) outputs[i] = output.data + i * cell.bytes;
        const auto encode = [&]() {
            Require(leo2_encode(codec.codec.get(), cell.bytes, inputs.data(), outputs.data(),
                scratch.data, scratch_bytes) == LEO2_SUCCESS, "encode");
        };
        const uint64_t input_hash = Hash(source.data, input_bytes);
        encode();
        const uint64_t output_hash = Hash(output.data, output_bytes);
        std::vector<long long> samples;
        if (measured)
        {
            for (unsigned i = 0; i < 4; ++i) encode();
            samples.reserve(21);
            for (unsigned i = 0; i < 21; ++i)
            {
                const auto start = std::chrono::steady_clock::now();
                encode();
                const auto end = std::chrono::steady_clock::now();
                samples.push_back(std::chrono::duration_cast<std::chrono::nanoseconds>(
                    end - start).count());
            }
        }
        Require(Hash(source.data, input_bytes) == input_hash &&
            Hash(output.data, output_bytes) == output_hash, "workload changed");
        if (argc == 5)
        {
            FILE* file = std::fopen(argv[4], "wbx");
            Require(file != NULL, "parity output exists or cannot be created");
            const bool ok = std::fwrite(output.data, 1, output_bytes, file) == output_bytes;
            const int close_status = std::fclose(file);
            Require(ok && close_status == 0, "parity write");
        }
        std::printf("{\"schema\":\"leopard-gfni-boundary/v1\",\"codec_commit\":\"%s\","
            "\"cell\":%u,\"k\":%u,\"r\":%u,\"bytes\":%zu,\"requested\":\"%s\","
            "\"execution_route\":\"%s\",\"scratch_bytes\":%zu,\"input_hash\":\"%016llx\","
            "\"output_hash\":\"%016llx\",\"samples_ns\":[", LEO_BOUNDARY_CODEC_COMMIT,
            index, cell.k, cell.r, cell.bytes, argv[3], backend == LEO2_BACKEND_GFNI ?
            "gfni" : "avx2", scratch_bytes, static_cast<unsigned long long>(input_hash),
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
#endif
