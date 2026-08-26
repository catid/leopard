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

namespace leopard {
void xor_mem_baseline(
    void* destination,
    const void* source,
    uint64_t byte_count);
}

namespace {

void require(bool condition, const char* message)
{
    if (!condition)
        throw std::runtime_error(message);
}

void test_next_power_of_two()
{
    // Keep the historical failure input runtime-visible so UBSan observes any
    // regression to calling LastNonzeroBit32(0), even in optimized builds.
    volatile unsigned one = 1;
    require(leopard::NextPow2(0) == 0,
        "NextPow2 zero sentinel");
    require(leopard::NextPow2(one) == 1,
        "NextPow2 one");
    require(leopard::NextPow2(2) == 2,
        "NextPow2 two");
    require(leopard::NextPow2(3) == 4,
        "NextPow2 rounds three");
    require(leopard::NextPow2(4) == 4,
        "NextPow2 four");
    require(leopard::NextPow2(5) == 8,
        "NextPow2 rounds five");
    require(leopard::NextPow2(0x7fffffffU) == 0x80000000U,
        "NextPow2 rounds to the largest representable power");
    require(leopard::NextPow2(0x80000000U) == 0x80000000U,
        "NextPow2 largest representable power");
    require(leopard::NextPow2(0x80000001U) == 0,
        "NextPow2 overflow sentinel");
    require(leopard::NextPow2(~0U) == 0,
        "NextPow2 maximum-input overflow sentinel");
}

void test_feature_classifier()
{
    using leopard::backend::ClassifyX86Features;
    require(!ClassifyX86Features(0, ~0U, ~0U, 0, ~0ULL).ssse3,
        "leaf-zero SSSE3 classification");
    require(!ClassifyX86Features(0, ~0U, ~0U, 0, ~0ULL).avx2,
        "leaf-zero AVX2 classification");
    require(!ClassifyX86Features(0, ~0U, ~0U, 0, ~0ULL).avx512,
        "leaf-zero AVX-512 classification");

    const uint32_t ssse3 = 0x00000200U;
    const uint32_t osxsave_avx = 0x18000000U;
    const uint32_t avx2 = 0x00000020U;
    const uint32_t avx512f = 0x00010000U;
    const uint32_t avx512bw = 0x40000000U;
    const uint32_t avx512vl = 0x80000000U;
    const uint32_t avx512 = avx512f | avx512bw | avx512vl;
    const uint32_t gfni = 0x00000100U;
    require(ClassifyX86Features(1, ssse3, 0, 0, 0).ssse3,
        "SSSE3 classification");
    require(!ClassifyX86Features(7, ssse3, avx2, 0, 6).avx2,
        "AVX2 without AVX/OSXSAVE");
    require(!ClassifyX86Features(7, ssse3 | osxsave_avx, avx2, 0, 2).avx2,
        "AVX2 without YMM state");
    require(!ClassifyX86Features(6, ssse3 | osxsave_avx, avx2, 0, 6).avx2,
        "AVX2 without leaf seven");
    const leopard::backend::X86Features complete = ClassifyX86Features(
        7, ssse3 | osxsave_avx, avx2, 0, 6);
    require(complete.ssse3 && complete.avx2,
        "complete AVX2 classification");
    require(!complete.avx512, "AVX2 misclassified as AVX-512");
    require(!complete.gfni, "AVX2 without the GFNI bit misclassified");
    require(!ClassifyX86Features(
            7, ssse3 | osxsave_avx, avx2 | avx512, 0, 6).avx512,
        "AVX-512 without opmask/ZMM state");
    require(!ClassifyX86Features(
            7, ssse3 | osxsave_avx, avx2 | avx512bw | avx512vl,
            0, 0xe6).avx512,
        "AVX-512 without AVX512F");
    require(!ClassifyX86Features(
            7, ssse3 | osxsave_avx, avx2 | avx512f | avx512vl,
            0, 0xe6).avx512,
        "AVX-512 without AVX512BW");
    require(!ClassifyX86Features(
            7, ssse3 | osxsave_avx, avx2 | avx512f | avx512bw,
            0, 0xe6).avx512,
        "AVX-512 without AVX512VL");
    const leopard::backend::X86Features complete_avx512 =
        ClassifyX86Features(
            7, ssse3 | osxsave_avx, avx2 | avx512, 0, 0xe6);
    require(complete_avx512.ssse3 && complete_avx512.avx2 &&
            complete_avx512.avx512,
        "complete AVX-512 classification");

    // GFNI is an AVX2-tier qualification: the affine kernels are VEX-encoded
    // and therefore require the same YMM state the AVX2 member requires, so
    // the leaf-seven ECX bit alone never qualifies the member.
    const leopard::backend::X86Features complete_gfni =
        ClassifyX86Features(
            7, ssse3 | osxsave_avx, avx2, gfni, 6);
    require(complete_gfni.ssse3 && complete_gfni.avx2 &&
            complete_gfni.gfni,
        "complete GFNI classification");
    require(!complete_gfni.avx512, "GFNI misclassified as AVX-512");
    require(!ClassifyX86Features(
            7, ssse3 | osxsave_avx, 0, gfni, 6).gfni,
        "GFNI without AVX2");
    require(!ClassifyX86Features(
            7, ssse3 | osxsave_avx, avx2, gfni, 2).gfni,
        "GFNI without YMM state");
}

void test_processor_classifier()
{
    static const uint32_t amd_ebx = 0x68747541U;
    static const uint32_t amd_edx = 0x69746e65U;
    static const uint32_t amd_ecx = 0x444d4163U;
    const uint32_t family_1a_model_44 =
        (0xfU << 8) | (0xbU << 20) | (4U << 4) | (4U << 16);
    const leopard::backend::X86ProcessorIdentity zen5 =
        leopard::backend::ClassifyX86Processor(
            amd_ebx, amd_edx, amd_ecx, family_1a_model_44);
    require(zen5.authentic_amd && zen5.family == 0x1aU &&
            zen5.model == 0x44U,
        "AMD extended family/model classification");
    require(leopard::backend::IsCalibratedAutoAVX512EncodeProcessor(zen5),
        "calibrated AMD processor classification");
    require(leopard::backend::IsCalibratedK1AVX2CopyProcessor(zen5),
        "calibrated K1 AVX2 copy processor classification");
    require(!leopard::backend::
            IsCalibratedK65R65B64AVX512GFNIProcessor(zen5),
        "Granite Ridge entered the Threadripper-only T128 selector");
    require(!leopard::backend::
            IsCalibratedK16R16B64AVX512GFNIProcessor(zen5),
        "Granite Ridge entered the Threadripper-only T16 selector");

    const uint32_t family_1a_model_08 =
        (0xfU << 8) | (0xbU << 20) | (8U << 4);
    const leopard::backend::X86ProcessorIdentity threadripper =
        leopard::backend::ClassifyX86Processor(
            amd_ebx, amd_edx, amd_ecx, family_1a_model_08);
    require(threadripper.authentic_amd && threadripper.family == 0x1aU &&
            threadripper.model == 0x08U,
        "Threadripper extended family/model classification");
    require(leopard::backend::
            IsCalibratedK65R65B64AVX512GFNIProcessor(threadripper),
        "calibrated Threadripper T128 processor classification");
    require(leopard::backend::
            IsCalibratedK16R16B64AVX512GFNIProcessor(threadripper),
        "calibrated Threadripper T16 processor classification");
    const uint32_t family_1a_model_07 =
        (0xfU << 8) | (0xbU << 20) | (7U << 4);
    const uint32_t family_1a_model_09 =
        (0xfU << 8) | (0xbU << 20) | (9U << 4);
    const leopard::backend::X86ProcessorIdentity adjacent_before =
        leopard::backend::ClassifyX86Processor(
            amd_ebx, amd_edx, amd_ecx, family_1a_model_07);
    const leopard::backend::X86ProcessorIdentity adjacent_after =
        leopard::backend::ClassifyX86Processor(
            amd_ebx, amd_edx, amd_ecx, family_1a_model_09);
    require(!leopard::backend::
                IsCalibratedK16R16B64AVX512GFNIProcessor(adjacent_before) &&
            !leopard::backend::
                IsCalibratedK16R16B64AVX512GFNIProcessor(adjacent_after),
        "adjacent AMD models entered the exact T16 selector");

    const leopard::backend::X86ProcessorIdentity wrong_vendor =
        leopard::backend::ClassifyX86Processor(
            amd_ebx ^ 1U, amd_edx, amd_ecx, family_1a_model_44);
    require(!wrong_vendor.authentic_amd &&
            wrong_vendor.family == 0x1aU && wrong_vendor.model == 0x44U,
        "vendor classification contaminated family/model");
    require(!leopard::backend::IsCalibratedAutoAVX512EncodeProcessor(
                wrong_vendor),
        "non-AMD processor entered the calibrated selector");
    require(!leopard::backend::IsCalibratedK1AVX2CopyProcessor(
                wrong_vendor),
        "non-AMD processor entered the calibrated K1 copy selector");
    const leopard::backend::X86ProcessorIdentity wrong_t128_vendor =
        leopard::backend::ClassifyX86Processor(
            amd_ebx ^ 1U, amd_edx, amd_ecx, family_1a_model_08);
    require(!leopard::backend::
            IsCalibratedK65R65B64AVX512GFNIProcessor(
                wrong_t128_vendor),
        "non-AMD processor entered the calibrated T128 selector");
    require(!leopard::backend::
            IsCalibratedK16R16B64AVX512GFNIProcessor(
                wrong_t128_vendor),
        "non-AMD processor entered the calibrated T16 selector");

    const uint32_t family_6_model_9e =
        (6U << 8) | (0xeU << 4) | (9U << 16);
    const leopard::backend::X86ProcessorIdentity family6 =
        leopard::backend::ClassifyX86Processor(
            0, 0, 0, family_6_model_9e);
    require(family6.family == 6U && family6.model == 0x9eU,
        "family-six extended model classification");
    require(!leopard::backend::IsCalibratedAutoAVX512EncodeProcessor(family6),
        "uncalibrated family/model entered the calibrated selector");
    require(!leopard::backend::IsCalibratedK1AVX2CopyProcessor(family6),
        "uncalibrated family/model entered the calibrated K1 copy selector");
    require(!leopard::backend::
            IsCalibratedK16R16B64AVX512GFNIProcessor(family6),
        "uncalibrated family/model entered the calibrated T16 selector");
}

#ifdef LEO_HAS_FF8
void test_avx2_six_output_large()
{
    const leopard::backend::Ops* ops =
        leopard::backend::GetQualifiedOps(LEO2_BACKEND_AVX2);
    if (!ops)
        return;
    require(ops->ff8_multiply_add_outputs != NULL,
        "qualified AVX2 backend omitted multi-output GF8 multiply-add");

    static const size_t byte_counts[] = {
        2047, 2048, 2049, 4096, 16384, 16385
    };
    static const uint16_t log_sets[][6] = {
        { 0, 1, 17, 29, 127, 254 },
        { 0, UINT16_MAX, 17, 29, UINT16_MAX, 254 }
    };
    static const size_t kStorageBytes = 16387;
    std::vector<uint8_t> source(kStorageBytes);
    std::vector<uint8_t> source_before(kStorageBytes);
    std::vector<std::vector<uint8_t> > outputs(
        6, std::vector<uint8_t>(kStorageBytes));
    std::vector<std::vector<uint8_t> > expected(
        6, std::vector<uint8_t>(kStorageBytes));
    void* output_pointers[6];

    for (size_t count_index = 0;
         count_index < sizeof(byte_counts) / sizeof(byte_counts[0]);
         ++count_index)
    {
        const size_t bytes = byte_counts[count_index];
        for (size_t log_set = 0;
             log_set < sizeof(log_sets) / sizeof(log_sets[0]);
             ++log_set)
        {
            for (size_t i = 0; i < kStorageBytes; ++i)
            {
                source[i] = source_before[i] = static_cast<uint8_t>(
                    i * 73U + count_index * 19U + log_set * 31U + 11U);
            }
            for (unsigned output = 0; output < 6; ++output)
            {
                output_pointers[output] = outputs[output].data() + 1;
                for (size_t i = 0; i < kStorageBytes; ++i)
                {
                    outputs[output][i] = expected[output][i] =
                        static_cast<uint8_t>(
                            i * (29U + output * 6U) +
                            count_index * 7U + log_set * 13U + output);
                }
                const uint16_t log = log_sets[log_set][output];
                if (log == UINT16_MAX)
                    continue;
                for (size_t i = 0; i < bytes; ++i)
                {
                    expected[output][i + 1] ^=
                        leopard::ff8::MultiplyLogElement(
                            source[i + 1], static_cast<uint8_t>(log));
                }
            }
            ops->ff8_multiply_add_outputs(
                output_pointers, source.data() + 1,
                log_sets[log_set], 6, bytes);
            require(source == source_before,
                "AVX2 six-output kernel modified its source");
            for (unsigned output = 0; output < 6; ++output)
            {
                require(outputs[output] == expected[output],
                    "AVX2 six-output large/tail result mismatch");
            }
        }
    }
}
#endif

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
    uint8_t original_x16[130];
    uint8_t original_y16[130];
    uint8_t accumulator_x16[130];
    uint8_t accumulator_y16[130];
    uint8_t expected_accumulator_x16[130];
    uint8_t expected_accumulator_y16[130];
    std::memcpy(butterfly_x16, source16, sizeof(source16));
    std::memcpy(butterfly_y16, product16, sizeof(product16));
    std::memcpy(original_x16, source16, sizeof(source16));
    std::memcpy(original_y16, product16, sizeof(product16));
    for (unsigned i = 0; i < sizeof(accumulator_x16); ++i)
    {
        accumulator_x16[i] = expected_accumulator_x16[i] =
            static_cast<uint8_t>(i * 61U + seed * 23U);
        accumulator_y16[i] = expected_accumulator_y16[i] =
            static_cast<uint8_t>(i * 103U + seed * 31U);
    }
    for (unsigned tile = 0; tile < 2; ++tile)
    {
        const unsigned offset = tile * 64U;
        for (unsigned lane = 0; lane < 32; ++lane)
        {
            uint16_t x_value = static_cast<uint16_t>(
                source16[offset + lane] |
                (static_cast<unsigned>(source16[offset + 32U + lane]) << 8));
            uint16_t y_value = static_cast<uint16_t>(
                product16[offset + lane] |
                (static_cast<unsigned>(product16[offset + 32U + lane]) << 8));
            y_value ^= x_value;
            x_value ^= leopard::ff16::MultiplyLogElement(y_value, log16);
            expected_accumulator_x16[offset + lane] ^=
                static_cast<uint8_t>(x_value);
            expected_accumulator_x16[offset + 32U + lane] ^=
                static_cast<uint8_t>(x_value >> 8);
            expected_accumulator_y16[offset + lane] ^=
                static_cast<uint8_t>(y_value);
            expected_accumulator_y16[offset + 32U + lane] ^=
                static_cast<uint8_t>(y_value >> 8);
        }
    }
    uint16_t tail_x = static_cast<uint16_t>(
        source16[128] | (static_cast<unsigned>(source16[129]) << 8));
    uint16_t tail_y = static_cast<uint16_t>(
        product16[128] | (static_cast<unsigned>(product16[129]) << 8));
    tail_y ^= tail_x;
    tail_x ^= leopard::ff16::MultiplyLogElement(tail_y, log16);
    expected_accumulator_x16[128] ^= static_cast<uint8_t>(tail_x);
    expected_accumulator_x16[129] ^= static_cast<uint8_t>(tail_x >> 8);
    expected_accumulator_y16[128] ^= static_cast<uint8_t>(tail_y);
    expected_accumulator_y16[129] ^= static_cast<uint8_t>(tail_y >> 8);

    ops.ff16_ifft_butterfly2_xor(
        original_x16, original_y16, accumulator_x16, accumulator_y16,
        log16, sizeof(original_x16));
    require(std::memcmp(accumulator_x16, expected_accumulator_x16,
                sizeof(accumulator_x16)) == 0 &&
            std::memcmp(accumulator_y16, expected_accumulator_y16,
                sizeof(accumulator_y16)) == 0,
        "concurrent GF16 accumulating butterfly mismatch");
    require(std::memcmp(original_x16, source16, sizeof(original_x16)) == 0 &&
            std::memcmp(original_y16, product16, sizeof(original_y16)) == 0,
        "concurrent GF16 accumulating butterfly modified an input");
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

void test_baseline_xor_exact_tails()
{
    static const size_t kGuardBytes = 79;
    static const size_t kMaximumBytes = 513;
    static const size_t kStorageBytes =
        kGuardBytes + kMaximumBytes + kGuardBytes;

    for (size_t bytes = 0; bytes <= kMaximumBytes; ++bytes)
    {
        std::vector<uint8_t> destination(kStorageBytes);
        std::vector<uint8_t> source(kStorageBytes);
        for (size_t i = 0; i < kStorageBytes; ++i)
        {
            destination[i] = static_cast<uint8_t>(
                i * 73U + bytes * 19U + 0x35U);
            source[i] = static_cast<uint8_t>(
                i * 29U + bytes * 47U + 0xa7U);
        }
        const std::vector<uint8_t> source_before = source;
        std::vector<uint8_t> expected = destination;
        for (size_t i = 0; i < bytes; ++i)
            expected[kGuardBytes + i] ^= source[kGuardBytes + i];

        leopard::xor_mem_baseline(
            destination.data() + kGuardBytes,
            source.data() + kGuardBytes,
            bytes);

        require(destination == expected,
            "scalar XOR changed a guard byte or mishandled an exact tail");
        require(source == source_before,
            "scalar XOR modified its read-only source");
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
    else if (std::strcmp(expected, "gfni") == 0 ||
             std::strcmp(expected, "avx2-gfni") == 0)
        require(backend == LEO2_BACKEND_GFNI, "forced GFNI selection");
    else
        throw std::runtime_error("invalid expected backend");
}

bool expected_forced_variant()
{
    const char* expected = std::getenv("LEO2_EXPECT_BACKEND");
    if (!expected || expected[0] == '\0')
        return false;
    return std::strcmp(expected, "scalar") == 0 ||
        std::strcmp(expected, "ssse3") == 0 ||
        std::strcmp(expected, "avx2") == 0 ||
        std::strcmp(expected, "avx512") == 0;
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
    // GFNI is the second explicit-only member and the same lower forced
    // variants cap it.  A forced AVX-512 variant does not cap it, and its
    // availability there still depends on the host, so it is not asserted.
    const leopard::backend::Ops* gfni =
        leopard::backend::GetQualifiedOps(LEO2_BACKEND_GFNI);
    if (forced_lower)
        require(gfni == NULL,
            "lower forced variant exposed explicit GFNI");
    else if (std::strcmp(expected, "gfni") == 0 ||
             std::strcmp(expected, "avx2-gfni") == 0)
        require(gfni != NULL,
            "forced GFNI variant did not qualify its selected table");
}

} // namespace

int main()
{
    try
    {
        test_next_power_of_two();
        test_feature_classifier();
        test_processor_classifier();
        require(leo_init() == Leopard_Success, "Leopard initialization");
        require(leopard::backend::StartupSelfTestPassed(),
            "backend startup known-answer tests");
#if defined(LEO2_HAVE_AVX512_GFNI_T128) && defined(LEO_HAS_FF8)
        {
            const leopard::backend::X86Features detected =
                leopard::backend::DetectX86Features();
            if (expected_forced_variant())
            {
                require(
                    leopard::backend::GetQualifiedAVX512GFNIT128() == NULL,
                    "forced variant published independent AVX-512/GFNI leaf");
                require(leopard::backend::
                        GetQualifiedAVX512GFNIT128Multiples64() == NULL,
                    "forced variant published larger AVX-512/GFNI leaf");
            }
            else
            {
                const bool expected_qualified = detected.avx512 &&
                    detected.gfni && leopard::backend::
                        IsCalibratedK65R65B64AVX512GFNIHost();
                require(
                    (leopard::backend::GetQualifiedAVX512GFNIT128() != NULL) ==
                        expected_qualified,
                    "AVX-512/GFNI T128 KAT publication disagrees with ISA gate");
                require((leopard::backend::
                            GetQualifiedAVX512GFNIT128Multiples64() != NULL) ==
                        expected_qualified,
                    "larger AVX-512/GFNI T128 KAT publication disagrees with ISA gate");
            }
        }
#else
        require(leopard::backend::GetQualifiedAVX512GFNIT128() == NULL,
            "uncompiled AVX-512/GFNI T128 kernel was published");
        require(leopard::backend::
                GetQualifiedAVX512GFNIT128Multiples64() == NULL,
            "uncompiled larger AVX-512/GFNI T128 kernel was published");
#endif
#if defined(LEO2_HAVE_AVX512_GFNI_T16) && defined(LEO_HAS_FF8)
        {
            const leopard::backend::X86Features detected =
                leopard::backend::DetectX86Features();
            if (expected_forced_variant())
            {
                require(leopard::backend::GetQualifiedAVX512GFNIT16() == NULL,
                    "forced variant published independent T16 AVX-512/GFNI leaf");
            }
            else
            {
                const bool expected_qualified = detected.avx512 &&
                    detected.gfni && leopard::backend::
                        IsCalibratedK16R16B64AVX512GFNIHost();
                require((leopard::backend::GetQualifiedAVX512GFNIT16() != NULL) ==
                        expected_qualified,
                    "AVX-512/GFNI T16 KAT publication disagrees with ISA gate");
            }
        }
#else
        require(leopard::backend::GetQualifiedAVX512GFNIT16() == NULL,
            "uncompiled AVX-512/GFNI T16 kernel was published");
#endif
        const char* baseline_xor_only =
            std::getenv("LEO2_TEST_BASELINE_XOR_ONLY");
        if (baseline_xor_only && std::strcmp(baseline_xor_only, "1") == 0)
        {
            test_baseline_xor_exact_tails();
            std::printf("Leopard2 baseline XOR exact-tail test passed\n");
            return 0;
        }
        const leopard::backend::Ops& ops = leopard::backend::GetOps();
        require(ops.kind == leopard::backend::SelectedBackend(),
            "fixed-ops introspection is not derived from selected ops");
        require(ops.kind != LEO2_BACKEND_AUTO && ops.name,
            "backend selection missing");
#ifdef LEO_HAS_FF8
        require((ops.kind == LEO2_BACKEND_AVX2 ||
                 ops.kind == LEO2_BACKEND_AVX512 ||
                 ops.kind == LEO2_BACKEND_GFNI) ==
                    (ops.ff8_weighted_ifft_butterfly4 != NULL),
            "weighted locator capability escaped the AVX2 GF8 backend");
#if defined(LEO2_GFNI_VARIANT) && !defined(LEO2_GFNI_MEMBER)
        const bool expected_ff8_walsh_locator = false;
#else
        const bool expected_ff8_walsh_locator =
            ops.kind == LEO2_BACKEND_AVX2;
#endif
        require(expected_ff8_walsh_locator ==
                    (ops.ff8_walsh_locator != NULL),
            "active locator capability escaped the pure AVX2 backend");
        const leopard::backend::Ops* scalar_ops =
            leopard::backend::GetQualifiedOps(LEO2_BACKEND_SCALAR);
        require(scalar_ops && !scalar_ops->ff8_weighted_ifft_butterfly4,
            "scalar backend exposed the AVX2-only weighted boundary");
        require(!scalar_ops->ff8_walsh_locator,
            "scalar backend exposed the AVX2-only active locator");
        const leopard::backend::Ops* ssse3_ops =
            leopard::backend::GetQualifiedOps(LEO2_BACKEND_SSSE3);
        if (ssse3_ops)
        {
            require(!ssse3_ops->ff8_weighted_ifft_butterfly4,
                "SSSE3 backend exposed the AVX2-only weighted boundary");
            require(!ssse3_ops->ff8_walsh_locator,
                "SSSE3 backend exposed the AVX2-only active locator");
        }
        const leopard::backend::Ops* avx2_ops =
            leopard::backend::GetQualifiedOps(LEO2_BACKEND_AVX2);
        if (avx2_ops)
        {
            require(avx2_ops->ff8_weighted_ifft_butterfly4 != NULL,
                "AVX2 backend omitted its qualified weighted boundary");
#if defined(LEO2_GFNI_VARIANT) && !defined(LEO2_GFNI_MEMBER)
            require(avx2_ops->ff8_walsh_locator == NULL,
                "in-place GFNI diagnostic exposed pure AVX2 locator");
#else
            require(avx2_ops->ff8_walsh_locator != NULL,
                "AVX2 backend omitted its qualified active locator");
#endif
        }
        const leopard::backend::Ops* avx512_ops =
            leopard::backend::GetQualifiedOps(LEO2_BACKEND_AVX512);
        if (avx512_ops)
        {
            require(avx512_ops->ff8_weighted_ifft_butterfly4 != NULL,
                "AVX-512 backend omitted its qualified weighted boundary");
            require(!avx512_ops->ff8_walsh_locator,
                "AVX-512 backend exposed the pure-AVX2 active locator");
        }
        const leopard::backend::Ops* gfni_ops =
            leopard::backend::GetQualifiedOps(LEO2_BACKEND_GFNI);
        if (gfni_ops)
        {
            require(gfni_ops->ff8_weighted_ifft_butterfly4 != NULL,
                "GFNI backend omitted its qualified weighted boundary");
            require(!gfni_ops->ff8_walsh_locator,
                "GFNI backend exposed the pure-AVX2 active locator");
        }
        test_avx2_six_output_large();
#else
        require(!ops.ff8_weighted_ifft_butterfly4,
            "FF8-disabled build exposed a weighted boundary");
        require(!ops.ff8_walsh_locator,
            "FF8-disabled build exposed an active locator callback");
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
        test_baseline_xor_exact_tails();
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
