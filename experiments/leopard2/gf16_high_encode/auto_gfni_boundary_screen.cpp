// leopard-79h.38.5.4.17.1: fixed dual-linked AUTO boundary qualification driver.
#ifdef LEO_AUTO_BOUNDARY_MAIN
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

#ifndef LEO_AUTO_BOUNDARY_CODEC_COMMIT
#error Supply the independently pinned codec source identity.
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
        void* raw = NULL;
        Require(posix_memalign(&raw, 64, bytes ? bytes : 64) == 0, "allocation");
        data = static_cast<uint8_t*>(raw);
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
struct Cell { unsigned k, r; size_t bytes; bool batch, explicit_avx2; };
const Cell kCells[] = {
    {1000, 200, 32768, false, false},
    {1000, 199, 65536, false, false},
    {1000, 200, 32768, true, false},
    {1000, 199, 65536, true, false},
    {1000, 200, 65536, false, false},
    {1000, 199, 32768, false, false},
    {1000, 200, 32768, false, true},
    {4096, 512, 4096, false, false}
};
}

int main(int argc, char** argv)
{
    try
    {
        Require((argc == 4 || argc == 5) && (!std::strcmp(argv[1], "--check") ||
            !std::strcmp(argv[1], "--measure")), "usage: --check|--measure cell[0..7] 0|1 [parity_file]");
        Require(std::strlen(argv[2]) == 1 && argv[2][0] >= '0' && argv[2][0] <= '7', "cell");
        Require(!std::strcmp(argv[3], "0") || !std::strcmp(argv[3], "1"), "boundary mode");
        const unsigned index = static_cast<unsigned>(argv[2][0] - '0');
        const Cell cell = kCells[index];
        const bool measured = !std::strcmp(argv[1], "--measure");
        Require(!measured || argc == 4, "no timed parity dump");
        const size_t input_bytes = static_cast<size_t>(cell.k) * cell.bytes;
        const size_t output_bytes = static_cast<size_t>(cell.r) * cell.bytes;
        Aligned source(input_bytes);
        std::vector<const void*> inputs(cell.k);
        uint32_t random = 20260906;
        for (size_t i = 0; i < input_bytes; ++i)
        {
            random ^= random << 13; random ^= random >> 17; random ^= random << 5;
            source.data[i] = static_cast<uint8_t>(random);
        }
        for (unsigned i = 0; i < cell.k; ++i) inputs[i] = source.data + i * cell.bytes;
#ifdef LEO_AUTO_BOUNDARY_MAIN
        Require(!std::strcmp(argv[3], "0"), "Leopard1 has no boundary mode");
        Require(leo_init() == 0 && leopard::CpuHasAVX2, "native Leopard1 initialization");
#ifndef LEO_TRY_AVX2
#error Native Leopard1 must contain AVX2.
#endif
        const int mode = -1;
        const char* implementation = "leopard1";
        const char* api = "leo_encode";
        const char* route = "native-avx2";
        const char* output_semantics = "first_r_work_buffers";
        const unsigned route_calls = 0;
        const unsigned count = leo_encode_work_count(cell.k, cell.r);
        Require(count >= cell.r, "Leopard1 scratch count");
        const size_t scratch_bytes = static_cast<size_t>(count) * cell.bytes;
        Aligned work(scratch_bytes);
        uint8_t* const parity = work.data;
        std::vector<void*> pointers(count);
        for (unsigned i = 0; i < count; ++i) pointers[i] = work.data + i * cell.bytes;
        const auto encode = [&]() {
            Require(leo_encode(cell.bytes, cell.k, cell.r, count, inputs.data(), pointers.data())
                == Leopard_Success, "Leopard1 encode");
        };
#else
        namespace diag = leopard2_internal;
        Require(!diag::AutoGF16GFNIBoundariesEnabledForDiagnostics(), "candidate must default off");
        const int mode = argv[3][0] - '0';
        Require(diag::SetAutoGF16GFNIBoundariesEnabledForDiagnostics(mode == 1), "boundary control");
        Require(diag::SetAutoGF16GFNIEncodeEnabledForDiagnostics(true), "arm untimed probe");
        leo2_context_options options = {};
        options.struct_size = sizeof(options);
        options.backend = cell.explicit_avx2 ? LEO2_BACKEND_AVX2 : LEO2_BACKEND_AUTO;
        options.thread_count = 1;
        leo2_context* raw_context = NULL;
        Require(leo2_context_create(&options, &raw_context) == LEO2_SUCCESS, "context");
        std::unique_ptr<leo2_context, decltype(&leo2_context_destroy)> context(raw_context, leo2_context_destroy);
        Require(leo2_context_backend(context.get()) == LEO2_BACKEND_AVX2, "context backend");
        leo2_codec* raw_codec = NULL;
        Require(leo2_codec_create(context.get(), cell.k, cell.r, LEO2_PROFILE_LEGACY_HIGH_V1,
            LEO2_FIELD_GF16, NULL, &raw_codec) == LEO2_SUCCESS, "codec");
        std::unique_ptr<leo2_codec, decltype(&leo2_codec_destroy)> codec(raw_codec, leo2_codec_destroy);
        const bool gfni = index == 4 || (index < 4 && mode == 1);
        Require(diag::AutoGF16GFNIEncodeSelectedForDiagnostics(codec.get(), cell.bytes) == gfni,
            "unexpected selected route");
        const unsigned route_calls = gfni ? 1 : 0;
        const char* implementation = "leopard2";
        const char* api = cell.batch ? "leo2_encode_batch_one_item" : "leo2_encode";
        const char* route = gfni ? "gfni" : "avx2";
        const char* output_semantics = "separate_output_buffers";
        size_t scratch_bytes = 0;
        Require(leo2_encode_scratch_size(codec.get(), cell.bytes, &scratch_bytes) == LEO2_SUCCESS,
            "scratch query");
        Aligned scratch(scratch_bytes), output(output_bytes);
        uint8_t* const parity = output.data;
        std::vector<void*> outputs(cell.r);
        for (unsigned i = 0; i < cell.r; ++i) outputs[i] = parity + i * cell.bytes;
        leo2_encode_batch_item item = {};
        item.shard_bytes = cell.bytes; item.original = inputs.data(); item.recovery = outputs.data();
        item.scratch = scratch.data; item.scratch_bytes = scratch_bytes;
        const auto encode = [&]() {
            Require((cell.batch ? leo2_encode_batch(codec.get(), &item, 1) :
                leo2_encode(codec.get(), cell.bytes, inputs.data(), outputs.data(), scratch.data,
                    scratch_bytes)) == LEO2_SUCCESS, "Leopard2 encode");
        };
#endif
        const uint64_t input_hash = Hash(source.data, input_bytes);
        encode();
        const uint64_t output_hash = Hash(parity, output_bytes);
#ifndef LEO_AUTO_BOUNDARY_MAIN
        Require(diag::AutoGF16GFNIEncodeCallCountForDiagnostics() == route_calls,
            "untimed public route count");
        Require(diag::FinishAutoGF16GFNIEncodeRouteProbeForDiagnostics() &&
            diag::AutoGF16GFNIEncodeModeForDiagnostics() == 1, "normalize probe before clocks");
#endif
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
                samples.push_back(std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count());
            }
        }
        Require(Hash(source.data, input_bytes) == input_hash && Hash(parity, output_bytes) == output_hash,
            "workload changed");
#ifndef LEO_AUTO_BOUNDARY_MAIN
        Require(diag::AutoGF16GFNIEncodeCallCountForDiagnostics() == route_calls &&
            diag::AutoGF16GFNIEncodeModeForDiagnostics() == 1 &&
            diag::AutoGF16GFNIBoundariesEnabledForDiagnostics() == (mode == 1) &&
            diag::AutoGF16GFNIEncodeSelectedForDiagnostics(codec.get(), cell.bytes) == gfni,
            "route or accounting changed");
#endif
        if (argc == 5)
        {
            FILE* file = std::fopen(argv[4], "wbx");
            Require(file != NULL, "parity file exists or cannot be created");
            const bool written = std::fwrite(parity, 1, output_bytes, file) == output_bytes;
            const int closed = std::fclose(file);
            Require(written && closed == 0, "parity write");
        }
        std::printf("{\"schema\":\"leopard-auto-gfni-boundary-screen/v1\",\"implementation\":\"%s\","
            "\"codec_commit\":\"%s\",\"cell\":%u,\"k\":%u,\"r\":%u,\"bytes\":%zu,"
            "\"boundary_mode\":%d,\"api\":\"%s\",\"execution_route\":\"%s\","
            "\"untimed_route_calls\":%u,\"scratch_bytes\":%zu,\"output_semantics\":\"%s\","
            "\"input_hash\":\"%016llx\",\"output_hash\":\"%016llx\",\"samples_ns\":[",
            implementation, LEO_AUTO_BOUNDARY_CODEC_COMMIT, index, cell.k, cell.r, cell.bytes,
            mode, api, route, route_calls, scratch_bytes, output_semantics,
            static_cast<unsigned long long>(input_hash), static_cast<unsigned long long>(output_hash));
        for (size_t i = 0; i < samples.size(); ++i) std::printf("%s%lld", i ? "," : "", samples[i]);
        std::puts("]}");
        return 0;
    }
    catch (const std::exception& error)
    {
        std::fprintf(stderr, "%s\n", error.what());
        return 1;
    }
}
