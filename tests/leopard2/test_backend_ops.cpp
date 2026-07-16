#include "Leopard2Backend.h"
#include "LeopardFF8.h"
#include "LeopardFF16.h"
#include "leopard.h"
#include "leopard2.h"

#include <atomic>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <stdexcept>
#include <thread>
#include <vector>

namespace {

void require(bool condition, const char* message)
{
    if (!condition)
        throw std::runtime_error(message);
}

void test_feature_classifier()
{
    using leopard::backend::ClassifyX86Features;
    require(!ClassifyX86Features(0, ~0U, ~0U, ~0ULL).ssse3,
        "leaf-zero SSSE3 classification");
    require(!ClassifyX86Features(0, ~0U, ~0U, ~0ULL).avx2,
        "leaf-zero AVX2 classification");

    const uint32_t ssse3 = 0x00000200U;
    const uint32_t osxsave_avx = 0x18000000U;
    const uint32_t avx2 = 0x00000020U;
    require(ClassifyX86Features(1, ssse3, 0, 0).ssse3,
        "SSSE3 classification");
    require(!ClassifyX86Features(7, ssse3, avx2, 6).avx2,
        "AVX2 without AVX/OSXSAVE");
    require(!ClassifyX86Features(7, ssse3 | osxsave_avx, avx2, 2).avx2,
        "AVX2 without YMM state");
    require(!ClassifyX86Features(6, ssse3 | osxsave_avx, avx2, 6).avx2,
        "AVX2 without leaf seven");
    const leopard::backend::X86Features complete = ClassifyX86Features(
        7, ssse3 | osxsave_avx, avx2, 6);
    require(complete.ssse3 && complete.avx2,
        "complete AVX2 classification");
}

void run_kernel_check(unsigned seed)
{
    const leopard::backend::Ops& ops = leopard::backend::GetOps();
    uint8_t source8[259];
    uint8_t product8[259];
    uint8_t expected8[259];
    for (unsigned i = 0; i < sizeof(source8); ++i)
        source8[i] = static_cast<uint8_t>(i * 37U + seed);
    const uint8_t log8 = static_cast<uint8_t>(seed * 29U);
    for (unsigned i = 0; i < sizeof(source8); ++i)
        expected8[i] = leopard::ff8::MultiplyLogElement(source8[i], log8);
    ops.ff8_multiply(product8, source8, log8, sizeof(source8));
    require(std::memcmp(product8, expected8, sizeof(product8)) == 0,
        "concurrent FF8 fixed multiply mismatch");

    uint8_t source16[130];
    uint8_t product16[130];
    uint8_t expected16[130];
    for (unsigned tile = 0; tile < 2; ++tile)
    {
        for (unsigned i = 0; i < 32; ++i)
        {
            const uint16_t value = static_cast<uint16_t>(
                (tile * 32U + i) * 23819U + seed);
            source16[tile * 64U + i] = static_cast<uint8_t>(value);
            source16[tile * 64U + 32U + i] =
                static_cast<uint8_t>(value >> 8);
        }
    }
    const uint16_t tail = static_cast<uint16_t>(seed * 40503U + 17U);
    source16[128] = static_cast<uint8_t>(tail);
    source16[129] = static_cast<uint8_t>(tail >> 8);
    const uint16_t log16 = static_cast<uint16_t>(seed * 7919U);
    for (unsigned tile = 0; tile < 2; ++tile)
    {
        for (unsigned i = 0; i < 32; ++i)
        {
            const uint16_t value = static_cast<uint16_t>(
                source16[tile * 64U + i] |
                (static_cast<unsigned>(source16[tile * 64U + 32U + i]) << 8));
            const uint16_t product =
                leopard::ff16::MultiplyLogElement(value, log16);
            expected16[tile * 64U + i] = static_cast<uint8_t>(product);
            expected16[tile * 64U + 32U + i] =
                static_cast<uint8_t>(product >> 8);
        }
    }
    const uint16_t tail_product =
        leopard::ff16::MultiplyLogElement(tail, log16);
    expected16[128] = static_cast<uint8_t>(tail_product);
    expected16[129] = static_cast<uint8_t>(tail_product >> 8);
    ops.ff16_multiply(product16, source16, log16, sizeof(source16));
    require(std::memcmp(product16, expected16, sizeof(product16)) == 0,
        "concurrent FF16 fixed multiply mismatch");
}

void test_concurrent_immutable_ops()
{
    std::atomic<unsigned> failures(0);
    std::vector<std::thread> threads;
    for (unsigned thread = 0; thread < 16; ++thread)
    {
        threads.push_back(std::thread([thread, &failures]() {
            try
            {
                for (unsigned iteration = 0; iteration < 64; ++iteration)
                    run_kernel_check(thread * 67U + iteration);
            }
            catch (...)
            {
                failures.fetch_add(1, std::memory_order_relaxed);
            }
        }));
    }
    for (size_t i = 0; i < threads.size(); ++i)
        threads[i].join();
    require(failures.load(std::memory_order_relaxed) == 0,
        "concurrent immutable backend execution failed");
}

void verify_expected_backend(leo2_backend backend)
{
    const char* expected = std::getenv("LEO2_EXPECT_BACKEND");
    if (!expected || expected[0] == '\0')
        return;
    if (std::strcmp(expected, "scalar") == 0)
        require(backend == LEO2_BACKEND_SCALAR, "forced scalar selection");
    else if (std::strcmp(expected, "ssse3") == 0)
        require(backend == LEO2_BACKEND_SSSE3, "forced SSSE3 selection");
    else if (std::strcmp(expected, "avx2") == 0)
        require(backend == LEO2_BACKEND_AVX2, "forced AVX2 selection");
    else
        throw std::runtime_error("invalid expected backend");
}

} // namespace

int main()
{
    try
    {
        test_feature_classifier();
        require(leo_init() == Leopard_Success, "Leopard initialization");
        require(leopard::backend::StartupSelfTestPassed(),
            "backend startup known-answer tests");
        const leopard::backend::Ops& ops = leopard::backend::GetOps();
        require(ops.kind == leopard::backend::SelectedBackend(),
            "backend introspection is not derived from selected ops");
        require(ops.kind != LEO2_BACKEND_AUTO && ops.name,
            "backend selection missing");

        leo2_context_options options;
        std::memset(&options, 0, sizeof(options));
        options.struct_size = sizeof(options);
        options.thread_count = 1;
        leo2_context* context = NULL;
        require(leo2_context_create(&options, &context) == LEO2_SUCCESS,
            "context creation");
        require(leo2_context_backend(context) == ops.kind,
            "public introspection differs from selected ops");
        verify_expected_backend(ops.kind);
        leo2_context_destroy(context);

        test_concurrent_immutable_ops();
        std::printf("Leopard2 backend ops passed: threads=16 "
            "iterations=1024 startup_kat=pass\n");
        return 0;
    }
    catch (const std::exception& error)
    {
        std::fprintf(stderr, "Leopard2 backend ops failed: %s\n", error.what());
        return 1;
    }
}
