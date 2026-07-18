#include "Leopard2Backend.h"
#include "LeopardCommon.h"
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
    require(!ClassifyX86Features(0, ~0U, ~0U, ~0ULL).avx512,
        "leaf-zero AVX-512 classification");

    const uint32_t ssse3 = 0x00000200U;
    const uint32_t osxsave_avx = 0x18000000U;
    const uint32_t avx2 = 0x00000020U;
    const uint32_t avx512f = 0x00010000U;
    const uint32_t avx512bw = 0x40000000U;
    const uint32_t avx512vl = 0x80000000U;
    const uint32_t avx512 = avx512f | avx512bw | avx512vl;
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
    require(!complete.avx512, "AVX2 misclassified as AVX-512");
    require(!ClassifyX86Features(
            7, ssse3 | osxsave_avx, avx2 | avx512, 6).avx512,
        "AVX-512 without opmask/ZMM state");
    require(!ClassifyX86Features(
            7, ssse3 | osxsave_avx, avx2 | avx512bw | avx512vl,
            0xe6).avx512,
        "AVX-512 without AVX512F");
    require(!ClassifyX86Features(
            7, ssse3 | osxsave_avx, avx2 | avx512f | avx512vl,
            0xe6).avx512,
        "AVX-512 without AVX512BW");
    require(!ClassifyX86Features(
            7, ssse3 | osxsave_avx, avx2 | avx512f | avx512bw,
            0xe6).avx512,
        "AVX-512 without AVX512VL");
    const leopard::backend::X86Features complete_avx512 =
        ClassifyX86Features(
            7, ssse3 | osxsave_avx, avx2 | avx512, 0xe6);
    require(complete_avx512.ssse3 && complete_avx512.avx2 &&
            complete_avx512.avx512,
        "complete AVX-512 classification");
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

    uint8_t butterfly_x8[259];
    uint8_t butterfly_y8[259];
    uint8_t original_x8[259];
    uint8_t original_y8[259];
    uint8_t accumulator_x8[259];
    uint8_t accumulator_y8[259];
    uint8_t expected_accumulator_x8[259];
    uint8_t expected_accumulator_y8[259];
    for (unsigned i = 0; i < sizeof(butterfly_x8); ++i)
    {
        butterfly_x8[i] = original_x8[i] = static_cast<uint8_t>(
            i * 43U + seed * 7U);
        butterfly_y8[i] = original_y8[i] = static_cast<uint8_t>(
            i * 83U + seed * 11U);
        accumulator_x8[i] = expected_accumulator_x8[i] =
            static_cast<uint8_t>(i * 97U + seed * 13U);
        accumulator_y8[i] = expected_accumulator_y8[i] =
            static_cast<uint8_t>(i * 17U + seed * 19U);
    }
    ops.ff8_ifft_butterfly2(
        butterfly_x8, butterfly_y8, log8, sizeof(butterfly_x8));
    for (unsigned i = 0; i < sizeof(butterfly_x8); ++i)
    {
        expected_accumulator_x8[i] ^= butterfly_x8[i];
        expected_accumulator_y8[i] ^= butterfly_y8[i];
    }
    ops.ff8_ifft_butterfly2_xor(
        original_x8, original_y8, accumulator_x8, accumulator_y8,
        log8, sizeof(original_x8));
    require(std::memcmp(accumulator_x8, expected_accumulator_x8,
                sizeof(accumulator_x8)) == 0 &&
            std::memcmp(accumulator_y8, expected_accumulator_y8,
                sizeof(accumulator_y8)) == 0,
        "concurrent GF8 accumulating butterfly mismatch");
    require(std::memcmp(original_x8, butterfly_x8,
                sizeof(original_x8)) != 0 ||
            std::memcmp(original_y8, butterfly_y8,
                sizeof(original_y8)) != 0,
        "concurrent GF8 butterfly fixture did not transform");
    ops.ff8_fft_butterfly2(
        butterfly_x8, butterfly_y8, log8, sizeof(butterfly_x8));
    require(std::memcmp(butterfly_x8, original_x8,
                sizeof(butterfly_x8)) == 0 &&
            std::memcmp(butterfly_y8, original_y8,
                sizeof(butterfly_y8)) == 0,
        "concurrent GF8 butterfly round trip mismatch");

    uint8_t out_x8[259];
    uint8_t out_y8[259];
    uint8_t expected_out_x8[259];
    uint8_t expected_out_y8[259];
    const uint8_t out_log8 = log8 == 255 ? 254 : log8;
    std::memcpy(expected_out_x8, original_x8, sizeof(original_x8));
    std::memcpy(expected_out_y8, original_y8, sizeof(original_y8));
    ops.ff8_fft_butterfly2(expected_out_x8, expected_out_y8,
        out_log8, sizeof(expected_out_x8));
    ops.ff8_fft_butterfly2_out(original_x8, original_y8, out_x8, out_y8,
        out_log8, sizeof(out_x8));
    require(std::memcmp(out_x8, expected_out_x8, sizeof(out_x8)) == 0 &&
            std::memcmp(out_y8, expected_out_y8, sizeof(out_y8)) == 0,
        "concurrent GF8 out-of-place butterfly mismatch");

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

    uint8_t butterfly_x16[130];
    uint8_t butterfly_y16[130];
    std::memcpy(butterfly_x16, source16, sizeof(source16));
    std::memcpy(butterfly_y16, product16, sizeof(product16));
    ops.ff16_ifft_butterfly2(
        butterfly_x16, butterfly_y16, log16, sizeof(butterfly_x16));
    ops.ff16_fft_butterfly2(
        butterfly_x16, butterfly_y16, log16, sizeof(butterfly_x16));
    require(std::memcmp(butterfly_x16, source16, sizeof(source16)) == 0 &&
            std::memcmp(butterfly_y16, product16, sizeof(product16)) == 0,
        "concurrent GF16 butterfly round trip mismatch");

    uint8_t out_x16[130];
    uint8_t out_y16[130];
    uint8_t expected_out_x16[130];
    uint8_t expected_out_y16[130];
    const uint16_t out_log16 = log16 == 65535 ? 65534 : log16;
    std::memcpy(expected_out_x16, source16, sizeof(source16));
    std::memcpy(expected_out_y16, product16, sizeof(product16));
    ops.ff16_fft_butterfly2(expected_out_x16, expected_out_y16,
        out_log16, sizeof(expected_out_x16));
    ops.ff16_fft_butterfly2_out(source16, product16, out_x16, out_y16,
        out_log16, sizeof(out_x16));
    require(std::memcmp(out_x16, expected_out_x16, sizeof(out_x16)) == 0 &&
            std::memcmp(out_y16, expected_out_y16, sizeof(out_y16)) == 0,
        "concurrent GF16 out-of-place butterfly mismatch");
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

void test_vector_xor_count_tail()
{
    static const unsigned kMaximumCount = 13;
    static const unsigned kBytes = 64;
    const leopard::backend::Ops& ops = leopard::backend::GetOps();

    for (unsigned threaded = 0; threaded < 2; ++threaded)
    {
        for (unsigned count = 1; count <= kMaximumCount; ++count)
        {
            std::vector<std::vector<uint8_t> > destinations(
                count, std::vector<uint8_t>(kBytes));
            std::vector<std::vector<uint8_t> > sources(
                count, std::vector<uint8_t>(kBytes));
            std::vector<std::vector<uint8_t> > expected(
                count, std::vector<uint8_t>(kBytes));
            std::vector<void*> destination_pointers(count);
            std::vector<void*> source_pointers(count);

            for (unsigned lane = 0; lane < count; ++lane)
            {
                for (unsigned byte = 0; byte < kBytes; ++byte)
                {
                    destinations[lane][byte] = static_cast<uint8_t>(
                        lane * 47U + byte * 13U + threaded);
                    sources[lane][byte] = static_cast<uint8_t>(
                        lane * 29U + byte * 71U + count);
                    expected[lane][byte] = static_cast<uint8_t>(
                        destinations[lane][byte] ^ sources[lane][byte]);
                }
                destination_pointers[lane] = destinations[lane].data();
                source_pointers[lane] = sources[lane].data();
            }

            if (threaded)
            {
                leopard::VectorXOR_Threads(
                    ops, kBytes, count,
                    destination_pointers.data(), source_pointers.data());
            }
            else
            {
                leopard::VectorXOR(
                    ops, kBytes, count,
                    destination_pointers.data(), source_pointers.data());
            }

            for (unsigned lane = 0; lane < count; ++lane)
            {
                require(destinations[lane] == expected[lane],
                    "VectorXOR count-tail mismatch");
            }
        }
    }
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
    else if (std::strcmp(expected, "avx512") == 0)
        require(backend == LEO2_BACKEND_AVX512,
            "forced AVX-512 selection");
    else
        throw std::runtime_error("invalid expected backend");
}

void verify_forced_avx512_qualification_cap()
{
    const char* expected = std::getenv("LEO2_EXPECT_BACKEND");
    if (!expected || expected[0] == '\0')
        return;
    const bool forced_lower = std::strcmp(expected, "scalar") == 0 ||
        std::strcmp(expected, "ssse3") == 0 ||
        std::strcmp(expected, "avx2") == 0;
    const leopard::backend::Ops* avx512 =
        leopard::backend::GetQualifiedOps(LEO2_BACKEND_AVX512);
    if (forced_lower)
        require(avx512 == NULL,
            "lower forced variant exposed explicit AVX-512");
    else if (std::strcmp(expected, "avx512") == 0)
        require(avx512 != NULL,
            "forced AVX-512 variant did not qualify its selected table");
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
            "fixed-ops introspection is not derived from selected ops");
        require(ops.kind != LEO2_BACKEND_AUTO && ops.name,
            "backend selection missing");
#ifdef LEO_HAS_FF8
        require((ops.kind == LEO2_BACKEND_AVX2 ||
                 ops.kind == LEO2_BACKEND_AVX512) ==
                    (ops.ff8_weighted_ifft_butterfly4 != NULL),
            "weighted locator capability escaped the AVX2 GF8 backend");
        const leopard::backend::Ops* scalar_ops =
            leopard::backend::GetQualifiedOps(LEO2_BACKEND_SCALAR);
        require(scalar_ops && !scalar_ops->ff8_weighted_ifft_butterfly4,
            "scalar backend exposed the AVX2-only weighted boundary");
        const leopard::backend::Ops* ssse3_ops =
            leopard::backend::GetQualifiedOps(LEO2_BACKEND_SSSE3);
        if (ssse3_ops)
            require(!ssse3_ops->ff8_weighted_ifft_butterfly4,
                "SSSE3 backend exposed the AVX2-only weighted boundary");
        const leopard::backend::Ops* avx2_ops =
            leopard::backend::GetQualifiedOps(LEO2_BACKEND_AVX2);
        if (avx2_ops)
            require(avx2_ops->ff8_weighted_ifft_butterfly4 != NULL,
                "AVX2 backend omitted its qualified weighted boundary");
        const leopard::backend::Ops* avx512_ops =
            leopard::backend::GetQualifiedOps(LEO2_BACKEND_AVX512);
        if (avx512_ops)
            require(avx512_ops->ff8_weighted_ifft_butterfly4 != NULL,
                "AVX-512 backend omitted its qualified weighted boundary");
#else
        require(!ops.ff8_weighted_ifft_butterfly4,
            "FF8-disabled build exposed a weighted boundary");
#endif

        leo2_context_options options;
        std::memset(&options, 0, sizeof(options));
        options.struct_size = sizeof(options);
        options.thread_count = 1;
        leo2_context* context = NULL;
        require(leo2_context_create(&options, &context) == LEO2_SUCCESS,
            "context creation");
        const leo2_backend execution = leopard::backend::ExecutionBackend();
        require(execution != LEO2_BACKEND_AUTO,
            "effective execution backend missing");
        require(leo2_context_backend(context) == execution,
            "public introspection differs from effective execution backend");
        verify_expected_backend(execution);
        verify_forced_avx512_qualification_cap();
        leo2_context_destroy(context);

        test_vector_xor_count_tail();
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
