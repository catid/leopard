// leopard-79h.38.5.4.18.1: native/pure-L1 versus explicit AVX2 L2.
// New attribution attempt; no codec source modifications or promotion claims.
#if defined(LEO_ISA_SCREEN_L1)
#include "leopard.h"
#include "LeopardCommon.h"
#else
#include "leopard2.h"
#include "Leopard2Direct.h"
#endif
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <memory>
#include <stdexcept>
#include <vector>
#include <omp.h>

#if !defined(LEO_ISA_PROFILE) || !defined(LEO_ISA_CODEC_COMMIT)
#error Bind the independent codec identity and build profile.
#endif

namespace {
void Require(bool ok, const char* message)
{
    if (!ok) throw std::runtime_error(message);
}

struct Aligned
{
    uint8_t* allocation;
    uint8_t* data;
    size_t bytes;
    explicit Aligned(size_t size) : allocation(NULL), data(NULL), bytes(size)
    {
        void* raw = NULL;
        Require(posix_memalign(&raw, 64, size + 128) == 0, "allocation");
        allocation = static_cast<uint8_t*>(raw);
        data = allocation + 64;
        std::memset(allocation, 0xa5, size + 128);
        std::memset(data, 0, size);
    }
    ~Aligned() { std::free(allocation); }
    Aligned(const Aligned&) = delete;
    Aligned& operator=(const Aligned&) = delete;
    void Check() const
    {
        for (unsigned i = 0; i < 64; ++i)
            Require(allocation[i] == 0xa5 && data[bytes + i] == 0xa5, "outer guard");
    }
};

uint64_t Hash(const uint8_t* data, size_t size)
{
    uint64_t value = UINT64_C(14695981039346656037);
    for (size_t i = 0; i < size; ++i)
        value = (value ^ data[i]) * UINT64_C(1099511628211);
    return value;
}

struct Cell { unsigned k, r; size_t bytes; };
const Cell cells[] = {{1000, 200, 65536}, {1000, 200, 32768},
                      {1000, 199, 65536}, {4096, 512, 4096}};
}

int main(int argc, char** argv)
{
    try
    {
        Require((argc == 3 || argc == 4) &&
                (!std::strcmp(argv[1], "--check") ||
                 !std::strcmp(argv[1], "--exercise") ||
                 !std::strcmp(argv[1], "--measure")),
                "usage: --check|--exercise|--measure cell[0..3] [untimed_parity_file]");
        Require(std::strlen(argv[2]) == 1 && argv[2][0] >= '0' && argv[2][0] <= '3',
                "invalid cell");
        const unsigned index = static_cast<unsigned>(argv[2][0] - '0');
        const Cell& cell = cells[index];
        const bool measured = !std::strcmp(argv[1], "--measure");
        const bool exercise = !std::strcmp(argv[1], "--exercise");
        Require(!measured || argc == 3, "no timed parity output");
        omp_set_dynamic(0);
        omp_set_num_threads(1);
        const size_t input_bytes = static_cast<size_t>(cell.k) * cell.bytes;
        const size_t output_bytes = static_cast<size_t>(cell.r) * cell.bytes;
        Aligned source(input_bytes);
        uint32_t random = 20260906;
        for (size_t i = 0; i < input_bytes; ++i)
        {
            random ^= random << 13;
            random ^= random >> 17;
            random ^= random << 5;
            source.data[i] = static_cast<uint8_t>(random);
        }
        std::vector<const void*> inputs(cell.k);
        for (unsigned i = 0; i < cell.k; ++i)
            inputs[i] = source.data + static_cast<size_t>(i) * cell.bytes;
        unsigned encode_calls = 0;
#if defined(LEO_ISA_SCREEN_L1)
        Require(leo_init() == 0 && leopard::CpuHasAVX2, "L1 AVX2 initialization");
        const unsigned count = leo_encode_work_count(cell.k, cell.r);
        Require(count >= cell.r, "L1 work count");
        const size_t scratch_bytes = static_cast<size_t>(count) * cell.bytes;
        Aligned scratch(scratch_bytes);
        uint8_t* const parity = scratch.data;
        const char* semantics = "first_r_work_buffers";
        std::vector<void*> outputs(count);
        for (unsigned i = 0; i < count; ++i)
            outputs[i] = scratch.data + static_cast<size_t>(i) * cell.bytes;
        const auto encode = [&]() {
            Require(leo_encode(cell.bytes, cell.k, cell.r, count,
                               inputs.data(), outputs.data()) == Leopard_Success, "L1 encode");
            ++encode_calls;
        };
#else
        leo2_context_options options = {};
        options.struct_size = sizeof(options);
        options.backend = LEO2_BACKEND_AVX2;
        options.thread_count = 1;
        leo2_context* raw_context = NULL;
        Require(leo2_context_create(&options, &raw_context) == LEO2_SUCCESS, "L2 context");
        std::unique_ptr<leo2_context, decltype(&leo2_context_destroy)>
            context(raw_context, leo2_context_destroy);
        Require(leo2_context_backend(context.get()) == LEO2_BACKEND_AVX2, "explicit AVX2");
        leo2_codec* raw_codec = NULL;
        Require(leo2_codec_create(context.get(), cell.k, cell.r,
                LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF16, NULL, &raw_codec) ==
                LEO2_SUCCESS, "L2 codec");
        std::unique_ptr<leo2_codec, decltype(&leo2_codec_destroy)>
            codec(raw_codec, leo2_codec_destroy);
        Require(!leopard2_internal::AutoGF16GFNIEncodeSelectedForDiagnostics(
                    codec.get(), cell.bytes), "explicit AVX2 switched to GFNI");
        size_t scratch_bytes = 0;
        Require(leo2_encode_scratch_size(codec.get(), cell.bytes, &scratch_bytes) ==
                LEO2_SUCCESS, "L2 scratch");
        Aligned scratch(scratch_bytes), output(output_bytes);
        uint8_t* const parity = output.data;
        const char* semantics = "separate_output_buffers";
        std::vector<void*> outputs(cell.r);
        for (unsigned i = 0; i < cell.r; ++i)
            outputs[i] = output.data + static_cast<size_t>(i) * cell.bytes;
        const auto encode = [&]() {
            Require(leo2_encode(codec.get(), cell.bytes, inputs.data(), outputs.data(),
                               scratch.data, scratch_bytes) == LEO2_SUCCESS, "L2 encode");
            ++encode_calls;
        };
#endif
        const uint64_t input_hash = Hash(source.data, input_bytes);
        encode();
        const uint64_t parity_hash = Hash(parity, output_bytes);
        long long samples[21] = {};
        if (measured || exercise)
        {
            for (unsigned i = 0; i < 4; ++i) encode();
            for (unsigned i = 0; i < 21; ++i)
            {
                if (measured)
                {
                    const auto start = std::chrono::steady_clock::now();
                    encode();
                    const auto end = std::chrono::steady_clock::now();
                    samples[i] = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
                }
                else encode();
            }
        }
        Require(encode_calls == (measured || exercise ? 26U : 1U), "encode call count");
        source.Check();
        scratch.Check();
        Require(Hash(source.data, input_bytes) == input_hash &&
                Hash(parity, output_bytes) == parity_hash, "workload changed");
#if !defined(LEO_ISA_SCREEN_L1)
        output.Check();
        Require(!leopard2_internal::AutoGF16GFNIEncodeSelectedForDiagnostics(
                    codec.get(), cell.bytes), "explicit AVX2 route changed");
#endif
        if (argc == 4)
        {
            FILE* file = std::fopen(argv[3], "wbx");
            Require(file != NULL, "parity file exists or cannot be created");
            const bool ok = std::fwrite(parity, 1, output_bytes, file) == output_bytes;
            const int status = std::fclose(file);
            Require(ok && status == 0, "parity write");
        }
        std::printf("{\"schema\":\"leopard-avx2-isa-screen/v1\",\"profile\":\"%s\","
                    "\"codec_commit\":\"%s\",\"cell\":%u,\"k\":%u,\"r\":%u,"
                    "\"bytes\":%zu,\"scratch_bytes\":%zu,\"output_semantics\":\"%s\","
                    "\"input_hash\":\"%016llx\",\"output_hash\":\"%016llx\","
                    "\"outer_guards\":true,\"encode_calls\":%u,\"samples_ns\":[",
                    LEO_ISA_PROFILE, LEO_ISA_CODEC_COMMIT, index, cell.k, cell.r,
                    cell.bytes, scratch_bytes, semantics,
                    static_cast<unsigned long long>(input_hash),
                    static_cast<unsigned long long>(parity_hash), encode_calls);
        if (measured)
            for (unsigned i = 0; i < 21; ++i)
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
