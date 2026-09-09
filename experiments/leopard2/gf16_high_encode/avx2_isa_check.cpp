// Clock-free exact-L1 ISA attribution checks, leopard-79h.38.5.4.18.
// This executable deliberately has no measurement option or clock reads.
#include "leopard.h"
#include "LeopardCommon.h"

#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <stdexcept>
#include <vector>
#include <omp.h>

#ifndef LEO_ATTRIBUTION_PROFILE
#error Supply the independently verified archive profile.
#endif
#define LEO_TEXT_INNER(x) #x
#define LEO_TEXT(x) LEO_TEXT_INNER(x)

namespace {
void Require(bool pass, const char* message)
{
    if (!pass) throw std::runtime_error(message);
}

struct Aligned
{
    uint8_t* allocation = NULL;
    uint8_t* data = NULL;
    size_t bytes;
    explicit Aligned(size_t count) : bytes(count)
    {
        void* raw = NULL;
        Require(posix_memalign(&raw, 64, count + 128) == 0, "allocation failed");
        allocation = static_cast<uint8_t*>(raw);
        data = allocation + 64;
        std::memset(allocation, 0xa5, count + 128);
        std::memset(data, 0, count);
    }
    ~Aligned() { std::free(allocation); }
    Aligned(const Aligned&) = delete;
    Aligned& operator=(const Aligned&) = delete;
    void Check() const
    {
        for (unsigned i = 0; i < 64; ++i)
            Require(allocation[i] == 0xa5 && data[bytes + i] == 0xa5,
                    "outer buffer guard changed");
    }
};

uint64_t Hash(const uint8_t* data, size_t bytes)
{
    uint64_t result = UINT64_C(14695981039346656037);
    for (size_t i = 0; i < bytes; ++i)
        result = (result ^ data[i]) * UINT64_C(1099511628211);
    return result;
}

struct Cell { unsigned k, r; size_t bytes; };
const Cell cells[] = {
    {1000, 200, 65536}, {1000, 200, 65536}, {1000, 200, 65536},
    {1000, 200, 32768}, {1000, 199, 65536}, {4096, 512, 4096},
    {100, 20, 64}, {12, 3, 4096}
};
}

int main(int argc, char** argv)
{
    try
    {
        Require(argc == 4 && std::strcmp(argv[1], "--check") == 0,
                "usage: avx2_isa_check --check cell[0..7] new_parity_file");
        Require(std::strlen(argv[2]) == 1 && argv[2][0] >= '0' && argv[2][0] <= '7',
                "invalid cell");
        const unsigned index = static_cast<unsigned>(argv[2][0] - '0');
        const Cell& cell = cells[index];
        omp_set_dynamic(0);
        omp_set_num_threads(1);
        Require(leo_init() == 0 && leopard::CpuHasAVX2, "AVX2 initialization failed");
        const size_t source_bytes = static_cast<size_t>(cell.k) * cell.bytes;
        const size_t parity_bytes = static_cast<size_t>(cell.r) * cell.bytes;
        Aligned source(source_bytes);
        uint32_t random = 20260906;
        for (size_t i = 0; i < source_bytes; ++i)
        {
            random ^= random << 13;
            random ^= random >> 17;
            random ^= random << 5;
            source.data[i] = static_cast<uint8_t>(random);
        }
        std::vector<const void*> originals(cell.k);
        for (unsigned i = 0; i < cell.k; ++i)
            originals[i] = source.data + static_cast<size_t>(i) * cell.bytes;
        const unsigned count = leo_encode_work_count(cell.k, cell.r);
        Require(count >= cell.r, "invalid work count");
        Aligned work(static_cast<size_t>(count) * cell.bytes);
        std::vector<void*> pointers(count);
        for (unsigned i = 0; i < count; ++i)
            pointers[i] = work.data + static_cast<size_t>(i) * cell.bytes;
        const uint64_t input_hash = Hash(source.data, source_bytes);
        Require(leo_encode(cell.bytes, cell.k, cell.r, count,
                          originals.data(), pointers.data()) == Leopard_Success,
                "encode failed");
        source.Check();
        work.Check();
        Require(Hash(source.data, source_bytes) == input_hash, "input changed");
        FILE* file = std::fopen(argv[3], "wbx");
        Require(file != NULL, "parity file exists or cannot be created");
        const bool written = std::fwrite(work.data, 1, parity_bytes, file) == parity_bytes;
        const int closed = std::fclose(file);
        Require(written && closed == 0, "parity dump failed");
        std::printf("{\"schema\":\"leopard-avx2-isa-check/v1\",\"timings\":false,"
                    "\"profile\":\"%s\",\"cell\":%u,\"k\":%u,\"r\":%u,"
                    "\"bytes\":%zu,\"parity_bytes\":%zu,\"input_hash\":\"%016llx\","
                    "\"output_hash\":\"%016llx\",\"outer_guards\":true}\n",
                    LEO_TEXT(LEO_ATTRIBUTION_PROFILE), index, cell.k, cell.r,
                    cell.bytes, parity_bytes,
                    static_cast<unsigned long long>(input_hash),
                    static_cast<unsigned long long>(Hash(work.data, parity_bytes)));
        return 0;
    }
    catch (const std::exception& error)
    {
        std::fprintf(stderr, "%s\n", error.what());
        return 1;
    }
}
