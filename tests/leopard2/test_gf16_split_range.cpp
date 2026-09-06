// Bounded deterministic coverage for leopard-79h.38.5.4.9.
// Retained range/encoder parity coverage after rejecting the cache-block trial.
#include "Leopard2Backend.h"
#include "LeopardFF16.h"
#include "leopard.h"

#include <cstdio>
#include <cstring>
#include <stdexcept>
#include <vector>

namespace {

void Require(bool value, const char* message)
{
    if (!value)
        throw std::runtime_error(message);
}

struct Rows
{
    // Deliberately unaligned, with an exact allocation end for ASan to check
    // overreads/writes.  There is no padded readable final SIMD vector.
    static const size_t kPrefix = 17;
    std::vector<std::vector<uint8_t> > storage;
    std::vector<void*> pointers;
    std::vector<const void*> sources;

    Rows(unsigned count, size_t bytes, uint32_t seed)
        : storage(count), pointers(count), sources(count)
    {
        for (unsigned row = 0; row < count; ++row)
        {
            storage[row].resize(kPrefix + bytes, 0xa5);
            for (size_t i = 0; i < bytes; ++i)
            {
                seed ^= seed << 13;
                seed ^= seed >> 17;
                seed ^= seed << 5;
                storage[row][kPrefix + i] = static_cast<uint8_t>(seed);
            }
            pointers[row] = storage[row].data() + kPrefix;
            sources[row] = pointers[row];
        }
    }

    void CheckPrefix() const
    {
        for (size_t row = 0; row < storage.size(); ++row)
            for (size_t i = 0; i < kPrefix; ++i)
                Require(storage[row][i] == 0xa5, "row prefix changed");
    }

};

// Whole-codec checks need many large rows.  Use slabs so ASan does not round
// hundreds of independently allocated 64-KiB-plus-prefix rows into larger
// allocator classes.  The range checks above retain exact per-row ends.
struct EncoderRows
{
    std::vector<uint8_t> storage;
    std::vector<void*> pointers;
    std::vector<const void*> sources;

    EncoderRows(unsigned count, size_t bytes, uint32_t seed)
        : storage(static_cast<size_t>(count) * bytes),
          pointers(count), sources(count)
    {
        for (size_t i = 0; i < storage.size(); ++i)
        {
            seed ^= seed << 13;
            seed ^= seed >> 17;
            seed ^= seed << 5;
            storage[i] = static_cast<uint8_t>(seed);
        }
        for (unsigned row = 0; row < count; ++row)
        {
            pointers[row] = storage.data() + static_cast<size_t>(row) * bytes;
            sources[row] = pointers[row];
        }
    }

    uint64_t Hash() const
    {
        uint64_t hash = UINT64_C(14695981039346656037);
        for (uint8_t value : storage)
            hash = (hash ^ value) * UINT64_C(1099511628211);
        return hash;
    }
};

void CheckRanges(const leopard::backend::Ops& candidate,
    const leopard::backend::Ops& scalar)
{
    const size_t lengths[] = {
        0, 2, 62, 64, 66, 128, 8190, 8192, 8194, 8254,
        16384, 16386, 32768, 32770, 65536, 65598
    };
    const unsigned distances[] = {1, 4, 16, 64};
    unsigned cases = 0;
    for (size_t length : lengths)
    for (unsigned distance : distances)
    for (unsigned zero_mask = 0; zero_mask < 8; ++zero_mask)
    for (unsigned inverse = 0; inverse < 2; ++inverse)
    for (unsigned prefer_fused = 0; prefer_fused < 2; ++prefer_fused)
    {
        const uint16_t logs[] = {
            static_cast<uint16_t>((zero_mask & 1) ? 65535 : 0),
            static_cast<uint16_t>((zero_mask & 2) ? 65535 : 32768),
            static_cast<uint16_t>((zero_mask & 4) ? 65535 : 65534)
        };
        Rows actual(distance * 4, length, 123 + zero_mask);
        Rows expected(distance * 4, length, 123 + zero_mask);
        const auto range = inverse ? candidate.ff16_ifft_butterfly4_range
                                   : candidate.ff16_fft_butterfly4_range;
        // Small and one-pair cases also exercise the independent scalar
        // table arithmetic.  Larger ranges use unchanged single-group calls
        // from this exact AVX-512 operation table, never the candidate range.
        const leopard::backend::Ops& oracle =
            (length <= 128 || distance == 1) ? scalar : candidate;
        const auto single = inverse ? oracle.ff16_ifft_butterfly4
                                    : oracle.ff16_fft_butterfly4;
        range(actual.pointers.data(), distance,
            logs[0], logs[1], logs[2], length, prefer_fused != 0);
        for (unsigned i = 0; i < distance; ++i)
            single(expected.pointers[i], expected.pointers[i + distance],
                expected.pointers[i + distance * 2],
                expected.pointers[i + distance * 3],
                logs[0], logs[1], logs[2], length);
        Require(actual.storage == expected.storage, "range differential");
        actual.CheckPrefix();
        ++cases;
    }
    std::printf("range differential: %u cases passed\n", cases);
}

template<bool Inverse>
void UnblockedRange(void* const* work, unsigned distance,
    uint16_t log01, uint16_t log23, uint16_t log02,
    uint64_t bytes, bool)
{
    const leopard::backend::Ops& ops =
        *leopard::backend::GetQualifiedOps(LEO2_BACKEND_AVX512);
    const auto single = Inverse ? ops.ff16_ifft_butterfly4
                                : ops.ff16_fft_butterfly4;
    for (unsigned i = 0; i < distance; ++i)
        single(work[i], work[i + distance], work[i + 2 * distance],
            work[i + 3 * distance], log01, log23, log02, bytes);
}

void CheckEncoders(const leopard::backend::Ops& candidate)
{
    leopard::backend::Ops control = candidate;
    // These callbacks are mandatory for a non-default operation table.
    // Replace them with the unchanged per-group route, not null pointers.
    control.ff16_ifft_butterfly4_range = UnblockedRange<true>;
    control.ff16_fft_butterfly4_range = UnblockedRange<false>;
    struct Shape { unsigned k, r, side; size_t bytes; };
    const Shape shapes[] = {
        {255, 129, 256, 8192}, {257, 200, 256, 8256},
        {1000, 199, 256, 32768}, {1000, 200, 256, 32768},
        {1000, 200, 256, 65536}, {4096, 512, 512, 4096}
    };
    for (const Shape& shape : shapes)
    {
        EncoderRows source(shape.k, shape.bytes, 20260906);
        EncoderRows actual(shape.side * 2, shape.bytes, 19);
        EncoderRows expected(shape.r, shape.bytes, 23);
        const uint64_t source_hash = source.Hash();
        leopard::ff16::ReedSolomonEncodeWithSourcePolicy(
            control, shape.bytes, 65536, shape.k, shape.r, shape.r,
            shape.side, source.sources.data(), actual.pointers.data(), NULL);
        for (unsigned row = 0; row < shape.r; ++row)
            std::memcpy(expected.pointers[row], actual.pointers[row],
                shape.bytes);
        Require(source.Hash() == source_hash, "control changed source");
        std::memset(actual.storage.data(), 0x5a, actual.storage.size());
        leopard::ff16::ReedSolomonEncodeWithSourcePolicy(
            candidate, shape.bytes, 65536, shape.k, shape.r, shape.r,
            shape.side, source.sources.data(), actual.pointers.data(), NULL);
        for (unsigned row = 0; row < shape.r; ++row)
            Require(std::memcmp(actual.pointers[row], expected.pointers[row],
                shape.bytes) == 0, "full encoder parity differential");
        Require(source.Hash() == source_hash, "encoder changed source");
        std::printf("encoder K%u/R%u/B%zu: parity matched\n",
            shape.k, shape.r, shape.bytes);
    }
}

} // namespace

int main(int argc, char** argv)
{
    std::setvbuf(stdout, NULL, _IONBF, 0);
    try
    {
        Require(argc == 2 && (std::strcmp(argv[1], "--ranges") == 0 ||
            std::strcmp(argv[1], "--encoder") == 0),
            "usage: test_gf16_split_range --ranges|--encoder");
        Require(leo_init() == Leopard_Success, "initialization failed");
        const leopard::backend::Ops* candidate =
            leopard::backend::GetQualifiedOps(LEO2_BACKEND_AVX512);
        if (!candidate)
        {
            std::puts("SKIP: AVX-512VL backend unavailable");
            return 77;
        }
        const leopard::backend::Ops* scalar =
            leopard::backend::GetQualifiedOps(LEO2_BACKEND_SCALAR);
        Require(scalar && candidate->ff16_ifft_butterfly4_range &&
            candidate->ff16_fft_butterfly4_range, "missing operation table");
        // Separate processes keep range-allocation sanitizer caches out of
        // the full encoder's live set, without reducing either test matrix.
        if (std::strcmp(argv[1], "--ranges") == 0)
            CheckRanges(*candidate, *scalar);
        else
            CheckEncoders(*candidate);
        return 0;
    }
    catch (const std::exception& error)
    {
        std::fprintf(stderr, "%s\n", error.what());
        return 1;
    }
}
