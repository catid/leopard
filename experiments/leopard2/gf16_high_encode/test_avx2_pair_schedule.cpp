// Untimed correctness for leopard-79h.38.5.4.18.2; separate from promotion.
#define LEO_BOUNDARY_GUARD_NO_MAIN 1
#include "test_gfni_boundary.cpp"
#include "Leopard2Backend.h"
#include "leopard.h"
#include <algorithm>
#include <exception>
#include <thread>

namespace {
void PairKernels()
{
    Require(leo_init() == Leopard_Success, "initialize");
    const auto* candidate = leopard::backend::GetQualifiedOps(LEO2_BACKEND_AVX2);
    const auto* scalar = leopard::backend::GetQualifiedOps(LEO2_BACKEND_SCALAR);
    Require(candidate && scalar, "AVX2 and scalar required, no skip");
    unsigned count = 0;
    // All ordinary multiplier logs, independently evaluated by scalar tables.
    // The all-ones zero-skew sentinel is a transform-caller specialization.
    Guard x(64, 1), y(64, 1), sx(64, 1), sy(64, 1);
    uint32_t random = 20260909;
    for (unsigned log = 0; log < 65535; ++log)
    {
        Fill(x.data, 64, random); Fill(y.data, 64, random);
        std::memcpy(sx.data, x.data, 64); std::memcpy(sy.data, y.data, 64);
        candidate->ff16_ifft_butterfly2(x.data, y.data, log, 64);
        scalar->ff16_ifft_butterfly2(sx.data, sy.data, log, 64);
        Require(!std::memcmp(x.data, sx.data, 64) &&
                !std::memcmp(y.data, sy.data, 64), "exhaustive inverse pair");
        x.Check(); y.Check(); sx.Check(); sy.Check();
        ++count;
    }
    const size_t lengths[] = {0, 2, 30, 32, 62, 64, 66, 126, 128, 130,
        8190, 8192, 16384, 32768, 32770, 65536, 65598};
    for (size_t bytes : lengths)
    for (size_t offset : {size_t(0), size_t(1), size_t(17)})
    for (unsigned log : {0U, 1U, 255U, 256U, 32768U, 65534U})
    for (unsigned inverse = 0; inverse < 2; ++inverse)
    {
        Guard a(bytes, offset), b(bytes, offset), c(bytes, offset), d(bytes, offset);
        Fill(a.data, bytes, random); Fill(b.data, bytes, random);
        std::memcpy(c.data, a.data, bytes); std::memcpy(d.data, b.data, bytes);
        const auto actual = inverse ? candidate->ff16_ifft_butterfly2
                                    : candidate->ff16_fft_butterfly2;
        const auto reference = inverse ? scalar->ff16_ifft_butterfly2
                                       : scalar->ff16_fft_butterfly2;
        actual(a.data, b.data, log, bytes);
        reference(c.data, d.data, log, bytes);
        Require(!std::memcmp(a.data, c.data, bytes) &&
                !std::memcmp(b.data, d.data, bytes), "pair boundary differential");
        a.Check(); b.Check(); c.Check(); d.Check();
        ++count;
    }
    std::printf("pair kernel cases: %u\n", count);
}

void SplitKernels()
{
    Require(leo_init() == Leopard_Success, "initialize");
    const auto* candidate = leopard::backend::GetQualifiedOps(LEO2_BACKEND_AVX2);
    const auto* scalar = leopard::backend::GetQualifiedOps(LEO2_BACKEND_SCALAR);
    Require(candidate && scalar, "AVX2 and scalar required, no skip");
    uint32_t random = 20260909;
    unsigned count = 0;
    for (size_t bytes : {size_t(0), size_t(2), size_t(62), size_t(64),
                        size_t(66), size_t(128), size_t(130), size_t(32768)})
    for (unsigned distance : {1U, 4U})
    for (unsigned zero = 0; zero < 8; ++zero)
    for (unsigned fused = 0; fused < 2; ++fused)
    {
        std::vector<std::unique_ptr<Guard> > actual, expected;
        std::vector<void*> pointers(4 * distance);
        for (unsigned i = 0; i < 4 * distance; ++i)
        {
            actual.emplace_back(new Guard(bytes, 17));
            expected.emplace_back(new Guard(bytes, 17));
            Fill(actual.back()->data, bytes, random);
            std::memcpy(expected.back()->data, actual.back()->data, bytes);
            pointers[i] = actual.back()->data;
        }
        const uint16_t a = zero & 1 ? 65535 : 0;
        const uint16_t b = zero & 2 ? 65535 : 32768;
        const uint16_t c = zero & 4 ? 65535 : 65534;
        candidate->ff16_ifft_butterfly4_range(pointers.data(), distance,
            a, b, c, bytes, fused != 0);
        for (unsigned i = 0; i < distance; ++i)
            scalar->ff16_ifft_butterfly4(expected[i]->data,
                expected[i + distance]->data, expected[i + 2 * distance]->data,
                expected[i + 3 * distance]->data, a, b, c, bytes);
        for (unsigned i = 0; i < 4 * distance; ++i)
        {
            Require(!std::memcmp(actual[i]->data, expected[i]->data, bytes),
                "split range differential");
            actual[i]->Check(); expected[i]->Check();
        }
        ++count;
    }
    std::printf("split range cases: %u\n", count);
}

void RoundTrip(unsigned cell)
{
    Require(cell < 3, "roundtrip cell");
    const unsigned k = cell == 0 ? 1000 : 17, r = cell == 0 ? 199 : 7;
    const size_t bytes = cell == 0 ? 32768 : cell == 1 ? 130 : 65;
    const auto field = cell == 2 ? LEO2_FIELD_GF8 : LEO2_FIELD_GF16;
    Codec candidate(k, r, field, cell == 0 ? LEO2_BACKEND_AUTO : LEO2_BACKEND_AVX2);
    Codec oracle(k, r, field, LEO2_BACKEND_GFNI);
    if (cell == 0)
        Require(!leopard2_internal::AutoGF16GFNIEncodeSelectedForDiagnostics(
            candidate.codec.get(), bytes), "AUTO AVX2 neighbor switched to GFNI");
    size_t decode_bytes = 0;
    Require(leo2_decode_scratch_size(candidate.codec.get(), bytes, &decode_bytes) ==
        LEO2_SUCCESS, "decode scratch");
    Guard scratch(std::max(candidate.Scratch(bytes), decode_bytes));
    Aligned inputs(static_cast<size_t>(k) * bytes), parity(static_cast<size_t>(r) * bytes);
    Aligned expected(static_cast<size_t>(r) * bytes);
    Guard restored(3 * bytes, 1);
    uint32_t random = 20260909;
    Fill(inputs.data, static_cast<size_t>(k) * bytes, random);
    const auto source_hash = Hash(inputs.data, static_cast<size_t>(k) * bytes);
    std::vector<const void*> sources(k), recovery(r);
    std::vector<void*> outputs(r), reference(r), recovered(k, NULL);
    std::vector<uint8_t> present(k, 1), repair_present(r, 1);
    for (unsigned i = 0; i < k; ++i) sources[i] = inputs.data + i * bytes;
    for (unsigned i = 0; i < r; ++i)
    {
        outputs[i] = parity.data + i * bytes;
        recovery[i] = outputs[i]; reference[i] = expected.data + i * bytes;
    }
    Require(leo2_encode(oracle.codec.get(), bytes, sources.data(), reference.data(),
        scratch.data, scratch.bytes) == LEO2_SUCCESS, "roundtrip oracle");
    Require(leo2_encode(candidate.codec.get(), bytes, sources.data(), outputs.data(),
        scratch.data, scratch.bytes) == LEO2_SUCCESS, "roundtrip candidate");
    Require(!std::memcmp(parity.data, expected.data, static_cast<size_t>(r) * bytes),
        "roundtrip full parity");
    for (unsigned i = 0; i < 3; ++i)
    {
        const unsigned row = i == 0 ? 0 : i == 1 ? k / 2 : k - 1;
        present[row] = 0; sources[row] = NULL;
        recovered[row] = restored.data + i * bytes;
    }
    Require(leo2_decode(candidate.codec.get(), bytes, present.data(), repair_present.data(),
        sources.data(), recovery.data(), recovered.data(), scratch.data, scratch.bytes) ==
        LEO2_SUCCESS, "roundtrip decode");
    for (unsigned i = 0; i < k; ++i)
        if (!present[i]) Require(!std::memcmp(recovered[i], inputs.data + i * bytes, bytes),
            "restored data mismatch");
    Require(source_hash == Hash(inputs.data, static_cast<size_t>(k) * bytes), "input changed");
    Require(!std::memcmp(parity.data, expected.data, static_cast<size_t>(r) * bytes),
        "decode changed parity");
    restored.Check(); scratch.Check();
    std::printf("roundtrip cell %u passed\n", cell);
}
}

int main(int argc, char** argv)
{
    try
    {
        Require(argc == 2, "one test selector required");
        if (!std::strcmp(argv[1], "--pairs")) PairKernels();
        else if (!std::strcmp(argv[1], "--split")) SplitKernels();
        else if (!std::strcmp(argv[1], "--roundtrip"))
            for (unsigned cell = 0; cell < 3; ++cell) RoundTrip(cell);
        else if (!std::strcmp(argv[1], "--concurrent"))
        {
            RoundTrip(1); // Finish lazy table qualification before shared reads.
            std::exception_ptr failures[4];
            std::vector<std::thread> threads;
            for (unsigned i = 0; i < 4; ++i)
                threads.emplace_back([&, i]() {
                    try { for (unsigned n = 0; n < 4; ++n) RoundTrip(1 + n % 2); }
                    catch (...) { failures[i] = std::current_exception(); }
                });
            for (auto& thread : threads) thread.join();
            for (auto failure : failures) if (failure) std::rethrow_exception(failure);
            std::puts("four-thread both-field roundtrips passed");
        }
        else
        {
            Require(std::strlen(argv[1]) == 1 && argv[1][0] >= '0' && argv[1][0] <= '7',
                "invalid selector");
            // The unchanged GFNI backend supplies the reference; partial-output
            // calls execute the actual candidate AVX2 backend, not the reverse.
            Check(argv[1][0] - '0', LEO2_BACKEND_AVX2, LEO2_BACKEND_GFNI);
        }
        return 0;
    }
    catch (const std::exception& error)
    {
        std::fprintf(stderr, "%s\n", error.what());
        return 1;
    }
}
