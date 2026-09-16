// Clock-free public-API qualification for leopard-79h.57.12.1.
// This is not a benchmark and intentionally has no --measure option.
#ifdef LEO_NATIVE_RELEASE_BASELINE
#include "leopard.h"
#else
#include "leopard2.h"
#endif

#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <memory>
#include <stdexcept>
#include <vector>
#include <omp.h>

#ifndef LEO_NATIVE_RELEASE_CODEC_COMMIT
#error Supply the source commit of the linked codec archive.
#endif

namespace {
void Require(bool condition, const char* message)
{
    if (!condition) throw std::runtime_error(message);
}

struct Cell { const char* id; unsigned k, r; size_t bytes; bool avx2; };
const Cell kCells[] = {
    {"copy", 1, 1, 4096, false},
    {"small", 16, 8, 64, false},
    {"gf8-high", 240, 16, 65536, false},
    {"gf8-balanced", 128, 128, 65536, false},
    {"gf16-inflation", 200, 50, 65536, false},
    {"gf16-gfni-region", 1000, 200, 65536, false},
    {"gf16-explicit-avx2", 1000, 200, 65536, true},
    {"gf16-large", 4096, 512, 4096, false}
};

struct Aligned
{
    void* data;
    size_t bytes;
    explicit Aligned(size_t count) : data(NULL), bytes(count)
    {
        if (count)
            Require(posix_memalign(&data, 64, count) == 0, "allocation failed");
    }
    ~Aligned() { std::free(data); }
    void Poison(unsigned char value) { if (bytes) std::memset(data, value, bytes); }
    Aligned(const Aligned&) = delete;
    Aligned& operator=(const Aligned&) = delete;
};

uint64_t Next(uint64_t& state)
{
    state ^= state << 13;
    state ^= state >> 7;
    state ^= state << 17;
    return state;
}

uint64_t Hash(const void* pointer, size_t bytes)
{
    const uint8_t* data = static_cast<const uint8_t*>(pointer);
    uint64_t hash = UINT64_C(14695981039346656037);
    for (size_t i = 0; i < bytes; ++i)
        hash = (hash ^ data[i]) * UINT64_C(1099511628211);
    return hash;
}

void Dump(const char* path, const void* data, size_t bytes)
{
    FILE* file = std::fopen(path, "wbx");
    Require(file != NULL, "parity output already exists or cannot be created");
    const bool written = std::fwrite(data, 1, bytes, file) == bytes;
    const int closed = std::fclose(file);
    Require(written && closed == 0, "parity write failed");
}
} // namespace

int main(int argc, char** argv)
{
    try
    {
        Require((argc == 3 || argc == 4) && std::strcmp(argv[1], "--check") == 0,
                "usage: encode_probe --check cell[0..7] [new_parity_file]");
        Require(std::strlen(argv[2]) == 1 && argv[2][0] >= '0' && argv[2][0] <= '7',
                "invalid cell");
        const unsigned index = static_cast<unsigned>(argv[2][0] - '0');
        const Cell& cell = kCells[index];
        omp_set_dynamic(0);
        omp_set_num_threads(1);
        Require(omp_get_max_threads() == 1, "OpenMP thread count changed");

        const size_t input_bytes = cell.k * cell.bytes;
        const size_t output_bytes = cell.r * cell.bytes;
        Aligned source(input_bytes);
        uint64_t random = UINT64_C(20260916);
        uint8_t* source_bytes = static_cast<uint8_t*>(source.data);
        for (size_t i = 0; i < input_bytes; ++i)
            source_bytes[i] = static_cast<uint8_t>(Next(random) >> 56);
        std::vector<const void*> inputs(cell.k);
        for (unsigned i = 0; i < cell.k; ++i)
            inputs[i] = source_bytes + i * cell.bytes;

#ifdef LEO_NATIVE_RELEASE_BASELINE
        Require(leo_init() == 0, "Leopard1 initialization failed");
        const char* implementation = "leopard1-native";
        const char* request = "native";
        const char* layout = "parity_is_first_r_work_rows";
        const int context_backend = -1; // No Leopard2 backend identifier applies.
        const unsigned work_count = leo_encode_work_count(cell.k, cell.r);
        Require(work_count >= cell.r, "Leopard1 work count");
        const size_t workspace_bytes = work_count * cell.bytes;
        const size_t separate_output_bytes = 0;
        Aligned scratch(workspace_bytes);
        void* const parity = scratch.data;
        std::vector<void*> work(work_count);
        for (unsigned i = 0; i < work_count; ++i)
            work[i] = static_cast<uint8_t*>(scratch.data) + i * cell.bytes;
        const auto encode = [&](unsigned char poison) {
            scratch.Poison(poison);
            Require(leo_encode(cell.bytes, cell.k, cell.r, work_count,
                inputs.data(), work.data()) == Leopard_Success, "Leopard1 encode failed");
        };
#else
        const char* implementation = "leopard2";
        const char* request = cell.avx2 ? "avx2" : "auto";
        const char* layout = "separate_parity_and_scratch";
        leo2_context_options options = {};
        options.struct_size = sizeof(options);
        options.backend = cell.avx2 ? LEO2_BACKEND_AVX2 : LEO2_BACKEND_AUTO;
        options.thread_count = 1;
        leo2_context* raw_context = NULL;
        Require(leo2_context_create(&options, &raw_context) == LEO2_SUCCESS,
                "context creation failed");
        std::unique_ptr<leo2_context, decltype(&leo2_context_destroy)>
            context(raw_context, leo2_context_destroy);
        const int context_backend = static_cast<int>(leo2_context_backend(context.get()));
        Require(leo2_context_thread_count(context.get()) == 1,
                "context thread count changed");
        Require(!cell.avx2 || context_backend == LEO2_BACKEND_AVX2,
                "explicit AVX2 request not preserved");
        leo2_codec* raw_codec = NULL;
        Require(leo2_codec_create(context.get(), cell.k, cell.r,
            LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_AUTO, NULL, &raw_codec) == LEO2_SUCCESS,
            "codec creation failed");
        std::unique_ptr<leo2_codec, decltype(&leo2_codec_destroy)>
            codec(raw_codec, leo2_codec_destroy);
        size_t workspace_bytes = 0;
        Require(leo2_encode_scratch_size(codec.get(), cell.bytes, &workspace_bytes)
                == LEO2_SUCCESS, "scratch query failed");
        const size_t separate_output_bytes = output_bytes;
        Aligned scratch(workspace_bytes), output(output_bytes);
        void* const parity = output.data;
        std::vector<void*> outputs(cell.r);
        for (unsigned i = 0; i < cell.r; ++i)
            outputs[i] = static_cast<uint8_t*>(output.data) + i * cell.bytes;
        const auto encode = [&](unsigned char poison) {
            scratch.Poison(poison);
            output.Poison(poison);
            Require(leo2_encode(codec.get(), cell.bytes, inputs.data(), outputs.data(),
                scratch.data, scratch.bytes) == LEO2_SUCCESS, "Leopard2 encode failed");
        };
#endif
        encode(0xa5);
        std::vector<uint8_t> first(output_bytes);
        std::memcpy(first.data(), parity, output_bytes);
        encode(0x5a);
        Require(std::memcmp(first.data(), parity, output_bytes) == 0,
                "repeated encode differs");
        random = UINT64_C(20260916);
        for (size_t i = 0; i < input_bytes; ++i)
            Require(source_bytes[i] == static_cast<uint8_t>(Next(random) >> 56),
                    "input changed");
        if (argc == 4) Dump(argv[3], parity, output_bytes);
        std::printf("{\"schema\":\"leopard-native-release-encode-check/v1\","
            "\"codec_commit\":\"%s\",\"implementation\":\"%s\","
            "\"cell\":%u,\"id\":\"%s\",\"k\":%u,\"r\":%u,\"bytes\":%zu,"
            "\"requested_backend\":\"%s\",\"context_backend\":%d,\"threads\":1,"
            "\"input_bytes\":%zu,\"parity_bytes\":%zu,"
            "\"workspace_bytes\":%zu,\"separate_output_bytes\":%zu,"
            "\"output_layout\":\"%s\",\"input_unchanged\":true,"
            "\"repeated_encode_equal\":true,\"public_encode_calls\":2,"
            "\"input_hash\":\"%016llx\",\"parity_hash\":\"%016llx\"}\n",
            LEO_NATIVE_RELEASE_CODEC_COMMIT, implementation, index, cell.id,
            cell.k, cell.r, cell.bytes, request, context_backend, input_bytes,
            output_bytes, workspace_bytes, separate_output_bytes, layout,
            static_cast<unsigned long long>(Hash(source.data, input_bytes)),
            static_cast<unsigned long long>(Hash(parity, output_bytes)));
        return 0;
    }
    catch (const std::exception& error)
    {
        std::fprintf(stderr, "%s\n", error.what());
        return 1;
    }
}
