// Encode-only diagnostic for leopard-79h.38.5.4.10; not v19 qualification.
// Compile separately against exact Leopard1 or current production Leopard2.
#ifdef LEO_CURRENT_SCREEN_MAIN
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

#ifndef LEO_CURRENT_SCREEN_SOURCE_COMMIT
#error Supply the codec source commit independently of the harness revision.
#endif

namespace {
void Require(bool value, const char* message)
{
    if (!value) throw std::runtime_error(message);
}

struct Aligned
{
    uint8_t* data;
    explicit Aligned(size_t bytes) : data(NULL)
    {
        void* allocation = NULL;
        Require(posix_memalign(&allocation, 64, bytes ? bytes : 64) == 0,
                "allocation failed");
        data = static_cast<uint8_t*>(allocation);
        std::memset(data, 0, bytes);
    }
    ~Aligned() { std::free(data); }
    Aligned(const Aligned&) = delete;
    Aligned& operator=(const Aligned&) = delete;
};

uint64_t Hash(const uint8_t* data, size_t bytes)
{
    uint64_t hash = UINT64_C(14695981039346656037);
    for (size_t i = 0; i < bytes; ++i)
        hash = (hash ^ data[i]) * UINT64_C(1099511628211);
    return hash;
}

enum Route { Auto, AVX2, AVX512 };
struct Cell { unsigned k, r; size_t bytes; Route requested; };
const Cell kCells[] = {
    {1000, 200, 65536, Auto},
    {1000, 200, 65536, AVX2},
    {1000, 200, 65536, AVX512},
    {1000, 200, 32768, Auto},
    {1000, 199, 65536, Auto},
    {4096, 512, 4096, Auto}
};
}

int main(int argc, char** argv)
{
    try
    {
        Require((argc == 3 || argc == 4) &&
                (std::strcmp(argv[1], "--check") == 0 ||
                 std::strcmp(argv[1], "--measure") == 0),
                "usage: current_route_screen --check|--measure cell[0..5] [parity_file]");
        Require(std::strlen(argv[2]) == 1 && argv[2][0] >= '0' &&
                argv[2][0] <= '5', "invalid cell");
        const unsigned index = static_cast<unsigned>(argv[2][0] - '0');
        const Cell& cell = kCells[index];
        const bool measured = std::strcmp(argv[1], "--measure") == 0;
        Require(!measured || argc == 3, "parity dump is untimed-only");
        const size_t input_bytes = static_cast<size_t>(cell.k) * cell.bytes;
        const size_t output_bytes = static_cast<size_t>(cell.r) * cell.bytes;
        Aligned source(input_bytes);
        std::vector<const void*> inputs(cell.k);
        uint32_t random = 20260906;
        for (size_t i = 0; i < input_bytes; ++i)
        {
            random ^= random << 13;
            random ^= random >> 17;
            random ^= random << 5;
            source.data[i] = static_cast<uint8_t>(random);
        }
        for (unsigned i = 0; i < cell.k; ++i)
            inputs[i] = source.data + static_cast<size_t>(i) * cell.bytes;

#ifdef LEO_CURRENT_SCREEN_MAIN
        Require(leo_init() == 0, "Leopard1 initialization failed");
#ifndef LEO_TRY_AVX2
#error The native Leopard1 profile must contain AVX2.
#endif
        Require(leopard::CpuHasAVX2, "Leopard1 AVX2 route unavailable");
        const char* implementation = "leopard1";
        const char* route = "native-avx2";
        const char* output_semantics = "first_r_work_buffers";
        const unsigned work_count = leo_encode_work_count(cell.k, cell.r);
        Require(work_count >= cell.r, "Leopard1 work count");
        const size_t scratch_bytes = static_cast<size_t>(work_count) * cell.bytes;
        Aligned work(scratch_bytes);
        uint8_t* const parity = work.data;
        std::vector<void*> work_pointers(work_count);
        for (unsigned i = 0; i < work_count; ++i)
            work_pointers[i] = work.data + static_cast<size_t>(i) * cell.bytes;
        const auto encode = [&]() {
            Require(leo_encode(cell.bytes, cell.k, cell.r, work_count,
                inputs.data(), work_pointers.data()) == Leopard_Success,
                "Leopard1 encode failed");
        };
#else
        leo2_context_options options = {};
        options.struct_size = sizeof(options);
        options.backend = cell.requested == Auto ? LEO2_BACKEND_AUTO :
            cell.requested == AVX2 ? LEO2_BACKEND_AVX2 : LEO2_BACKEND_AVX512;
        options.thread_count = 1;
        leo2_context* raw_context = NULL;
        Require(leo2_context_create(&options, &raw_context) == LEO2_SUCCESS,
                "context creation failed");
        std::unique_ptr<leo2_context, decltype(&leo2_context_destroy)>
            context(raw_context, leo2_context_destroy);
        leo2_codec* raw_codec = NULL;
        Require(leo2_codec_create(context.get(), cell.k, cell.r,
            LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF16, NULL, &raw_codec)
                == LEO2_SUCCESS, "codec creation failed");
        std::unique_ptr<leo2_codec, decltype(&leo2_codec_destroy)>
            codec(raw_codec, leo2_codec_destroy);
        const bool gfni = leopard2_internal::
            AutoGF16GFNIEncodeSelectedForDiagnostics(codec.get(), cell.bytes);
        Require(gfni == (index == 0), "unexpected AUTO GFNI route");
        Require(leo2_context_backend(context.get()) ==
            (cell.requested == AVX512 ? LEO2_BACKEND_AVX512 : LEO2_BACKEND_AVX2),
            "unexpected context backend");
        const char* implementation = "leopard2";
        const char* route = gfni ? "gfni" :
            cell.requested == AVX512 ? "avx512" : "avx2";
        const char* output_semantics = "separate_output_buffers";
        size_t scratch_bytes = 0;
        Require(leo2_encode_scratch_size(codec.get(), cell.bytes, &scratch_bytes)
                == LEO2_SUCCESS, "scratch query failed");
        Aligned scratch(scratch_bytes), output(output_bytes);
        uint8_t* const parity = output.data;
        std::vector<void*> outputs(cell.r);
        for (unsigned i = 0; i < cell.r; ++i)
            outputs[i] = parity + static_cast<size_t>(i) * cell.bytes;
        const auto encode = [&]() {
            Require(leo2_encode(codec.get(), cell.bytes, inputs.data(),
                outputs.data(), scratch.data, scratch_bytes) == LEO2_SUCCESS,
                "Leopard2 encode failed");
        };
#endif
        const uint64_t input_hash = Hash(source.data, input_bytes);
        encode();
        const uint64_t output_hash = Hash(parity, output_bytes);
        std::vector<long long> samples;
        // One public full-output encode per sample. No setup/allocation/hash,
        // extra Leopard1 parity copy, batching, or amortized plan reuse clock.
        // --check never reads a clock and can dump every parity byte.
        if (measured)
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
        Require(Hash(source.data, input_bytes) == input_hash &&
                Hash(parity, output_bytes) == output_hash, "workload changed");
#ifndef LEO_CURRENT_SCREEN_MAIN
        Require(leopard2_internal::AutoGF16GFNIEncodeSelectedForDiagnostics(
            codec.get(), cell.bytes) == gfni, "AUTO route changed");
#endif
        if (argc == 4)
        {
            FILE* file = std::fopen(argv[3], "wbx");
            Require(file != NULL, "parity file exists or cannot be created");
            const bool written = std::fwrite(parity, 1, output_bytes, file) == output_bytes;
            const int closed = std::fclose(file);
            Require(written && closed == 0, "parity dump failed");
        }
        std::printf("{\"schema\":\"leopard-gf16-current-route-screen/v1\","
            "\"implementation\":\"%s\",\"codec_commit\":\"%s\","
            "\"cell\":%u,\"k\":%u,\"r\":%u,\"bytes\":%zu,"
            "\"execution_route\":\"%s\",\"scratch_bytes\":%zu,"
            "\"output_semantics\":\"%s\",\"input_hash\":\"%016llx\","
            "\"output_hash\":\"%016llx\",\"samples_ns\":[", implementation,
            LEO_CURRENT_SCREEN_SOURCE_COMMIT, index, cell.k, cell.r, cell.bytes,
            route, scratch_bytes, output_semantics,
            static_cast<unsigned long long>(input_hash),
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
