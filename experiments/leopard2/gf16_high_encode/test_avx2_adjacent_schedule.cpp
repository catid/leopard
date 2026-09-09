// Untimed semantic qualification; leopard-79h.38.5.4.18.3.
#define main ExistingPairCheckMain
#include "test_avx2_pair_schedule.cpp"
#undef main

namespace {
void AdjacentPairs()
{
    Require(leo_init() == Leopard_Success, "initialize");
    const auto* candidate = leopard::backend::GetQualifiedOps(LEO2_BACKEND_AVX2);
    const auto* scalar = leopard::backend::GetQualifiedOps(LEO2_BACKEND_SCALAR);
    Require(candidate && scalar, "AVX2 and scalar required, no skip");
    Guard x(64, 1), y(64, 17), sx(64, 1), sy(64, 17);
    Guard u(64, 1), v(64, 17), su(64, 1), sv(64, 17);
    uint32_t random = 20260909;
    unsigned forward = 0, accumulating = 0;
    for (unsigned log = 0; log < 65535; ++log)
    {
        Fill(x.data, 64, random); Fill(y.data, 64, random);
        std::memcpy(sx.data, x.data, 64); std::memcpy(sy.data, y.data, 64);
        candidate->ff16_fft_butterfly2(x.data, y.data, log, 64);
        scalar->ff16_fft_butterfly2(sx.data, sy.data, log, 64);
        Require(!std::memcmp(x.data, sx.data, 64) && !std::memcmp(y.data, sy.data, 64),
                "exhaustive forward pair");
        const auto x_hash = Hash(x.data, 64), y_hash = Hash(y.data, 64);
        Fill(u.data, 64, random); Fill(v.data, 64, random);
        std::memcpy(su.data, u.data, 64); std::memcpy(sv.data, v.data, 64);
        candidate->ff16_ifft_butterfly2_xor(x.data, y.data, u.data, v.data, log, 64);
        scalar->ff16_ifft_butterfly2_xor(sx.data, sy.data, su.data, sv.data, log, 64);
        Require(!std::memcmp(u.data, su.data, 64) && !std::memcmp(v.data, sv.data, 64),
                "exhaustive accumulating pair");
        Require(Hash(x.data, 64) == x_hash && Hash(y.data, 64) == y_hash &&
                !std::memcmp(x.data, sx.data, 64) && !std::memcmp(y.data, sy.data, 64),
                "accumulating inputs changed");
        for (Guard* guard : {&x, &y, &sx, &sy, &u, &v, &su, &sv}) guard->Check();
        ++forward; ++accumulating;
    }
    unsigned boundaries = 0;
    for (size_t bytes : {size_t(0), size_t(2), size_t(30), size_t(32), size_t(62), size_t(64),
                        size_t(66), size_t(126), size_t(128), size_t(130), size_t(8190),
                        size_t(8192), size_t(16384), size_t(32768), size_t(32770),
                        size_t(65536), size_t(65598)})
    for (size_t offset : {size_t(0), size_t(1), size_t(17)})
    for (unsigned log : {0U, 1U, 255U, 256U, 32768U, 65534U})
    {
        Guard a(bytes, offset), b(bytes, offset), c(bytes, offset), d(bytes, offset);
        Guard sa(bytes, offset), sb(bytes, offset), sc(bytes, offset), sd(bytes, offset);
        Fill(a.data, bytes, random); Fill(b.data, bytes, random);
        Fill(c.data, bytes, random); Fill(d.data, bytes, random);
        std::memcpy(sa.data, a.data, bytes); std::memcpy(sb.data, b.data, bytes);
        std::memcpy(sc.data, c.data, bytes); std::memcpy(sd.data, d.data, bytes);
        std::vector<uint8_t> original_c(c.data, c.data + bytes), original_d(d.data, d.data + bytes);
        for (unsigned iteration = 0; iteration < 3; ++iteration)
        {
            candidate->ff16_ifft_butterfly2_xor(a.data, b.data, c.data, d.data, log, bytes);
            scalar->ff16_ifft_butterfly2_xor(sa.data, sb.data, sc.data, sd.data, log, bytes);
            Require(!std::memcmp(c.data, sc.data, bytes) && !std::memcmp(d.data, sd.data, bytes),
                    "accumulating boundary differential");
            Require(!std::memcmp(a.data, sa.data, bytes) && !std::memcmp(b.data, sb.data, bytes),
                    "accumulating boundary input mutation");
            if (iteration == 1 && bytes)
                Require(!std::memcmp(c.data, original_c.data(), bytes) &&
                        !std::memcmp(d.data, original_d.data(), bytes), "XOR accumulation involution");
            for (Guard* guard : {&a, &b, &c, &d, &sa, &sb, &sc, &sd}) guard->Check();
            ++boundaries;
        }
    }
    std::printf("adjacent pairs: forward=%u accumulating=%u boundary_accumulations=%u\n",
                forward, accumulating, boundaries);
}

void ForwardRanges()
{
    Require(leo_init() == Leopard_Success, "initialize");
    const auto* candidate = leopard::backend::GetQualifiedOps(LEO2_BACKEND_AVX2);
    const auto* scalar = leopard::backend::GetQualifiedOps(LEO2_BACKEND_SCALAR);
    Require(candidate && scalar, "AVX2 and scalar required, no skip");
    uint32_t random = 20260909;
    unsigned count = 0;
    for (size_t bytes : {size_t(0), size_t(2), size_t(62), size_t(64), size_t(66),
                        size_t(128), size_t(130), size_t(32768)})
    for (unsigned distance : {1U, 4U, 16U})
    for (unsigned zero = 0; zero < 8; ++zero)
    for (unsigned fused = 0; fused < 2; ++fused)
    {
        std::vector<std::unique_ptr<Guard> > actual, expected;
        std::vector<void*> pointers(4 * distance);
        for (unsigned i = 0; i < 4 * distance; ++i)
        {
            actual.emplace_back(new Guard(bytes, 17)); expected.emplace_back(new Guard(bytes, 17));
            Fill(actual.back()->data, bytes, random);
            std::memcpy(expected.back()->data, actual.back()->data, bytes);
            pointers[i] = actual.back()->data;
        }
        const uint16_t a = zero & 1 ? 65535 : 0, b = zero & 2 ? 65535 : 32768,
                       c = zero & 4 ? 65535 : 65534;
        candidate->ff16_fft_butterfly4_range(pointers.data(), distance, a, b, c, bytes, fused != 0);
        for (unsigned i = 0; i < distance; ++i)
            scalar->ff16_fft_butterfly4(expected[i]->data, expected[i + distance]->data,
                expected[i + 2 * distance]->data, expected[i + 3 * distance]->data, a, b, c, bytes);
        for (unsigned i = 0; i < 4 * distance; ++i)
        {
            Require(!std::memcmp(actual[i]->data, expected[i]->data, bytes), "forward range oracle");
            actual[i]->Check(); expected[i]->Check();
        }
        ++count;
    }
    std::printf("forward range cases: %u\n", count);
}
}

int main(int argc, char** argv)
{
    try
    {
        if (argc == 2 && !std::strcmp(argv[1], "--adjacent")) AdjacentPairs();
        else if (argc == 2 && !std::strcmp(argv[1], "--forward-ranges")) ForwardRanges();
        else return ExistingPairCheckMain(argc, argv);
        return 0;
    }
    catch (const std::exception& error)
    {
        std::fprintf(stderr, "%s\n", error.what());
        return 1;
    }
}
