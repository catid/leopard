/*
    Copyright (c) 2017 Christopher A. Taylor.  All rights reserved.
    See Leopard2Backend.h for the full BSD license notice.
*/

#include "Leopard2Backend.h"
#include "Leopard2Direct.h"
#include "leopard.h"
#include "leopard2.h"

#include <atomic>
#include <cstdio>
#include <cstring>
#include <stdexcept>
#include <thread>
#include <vector>

#ifndef LEO2_ENABLE_TEST_HOOKS
#error "Backend failure tests require LEO2_ENABLE_TEST_HOOKS"
#endif

namespace {

void require(bool condition, const char* message)
{
    if (!condition)
        throw std::runtime_error(message);
}

leo2_backend parse_backend(const char* name)
{
    if (std::strcmp(name, "scalar") == 0)
        return LEO2_BACKEND_SCALAR;
    if (std::strcmp(name, "ssse3") == 0)
        return LEO2_BACKEND_SSSE3;
    if (std::strcmp(name, "avx2") == 0)
        return LEO2_BACKEND_AVX2;
    if (std::strcmp(name, "avx512") == 0)
        return LEO2_BACKEND_AVX512;
    if (std::strcmp(name, "gfni") == 0 ||
        std::strcmp(name, "avx2-gfni") == 0)
        return LEO2_BACKEND_GFNI;
    throw std::runtime_error("unknown backend argument");
}

leopard::backend::TestSetupFault fault_for(
    leo2_backend backend,
    const char* stage)
{
    using namespace leopard::backend;
    /*
        Every backend is matched explicitly.  An unmatched backend must throw
        rather than fall through to another member's fault: silently injecting
        an AVX-512 fault for an unrelated backend would report a passing test
        for a member that was never exercised.
    */
    if (std::strcmp(stage, "ff8-allocation") == 0)
    {
        if (backend == LEO2_BACKEND_SCALAR)
            return TestSetupFaultScalarFF8Allocation;
        if (backend == LEO2_BACKEND_SSSE3)
            return TestSetupFaultSSSE3FF8Allocation;
        if (backend == LEO2_BACKEND_AVX2)
            return TestSetupFaultAVX2FF8Allocation;
        if (backend == LEO2_BACKEND_AVX512)
            return TestSetupFaultAVX512FF8Allocation;
        if (backend == LEO2_BACKEND_GFNI)
            return TestSetupFaultGFNIFF8Allocation;
        throw std::runtime_error("unknown failure-backend argument");
    }
    if (std::strcmp(stage, "ff16-allocation") == 0)
    {
        if (backend == LEO2_BACKEND_SCALAR)
            return TestSetupFaultScalarFF16Allocation;
        if (backend == LEO2_BACKEND_SSSE3)
            return TestSetupFaultSSSE3FF16Allocation;
        if (backend == LEO2_BACKEND_AVX2)
            return TestSetupFaultAVX2FF16Allocation;
        if (backend == LEO2_BACKEND_AVX512)
            return TestSetupFaultAVX512FF16Allocation;
        if (backend == LEO2_BACKEND_GFNI)
            return TestSetupFaultGFNIFF16Allocation;
        throw std::runtime_error("unknown failure-backend argument");
    }
    if (std::strcmp(stage, "kat") == 0)
    {
        if (backend == LEO2_BACKEND_SCALAR)
            return TestSetupFaultScalarKAT;
        if (backend == LEO2_BACKEND_SSSE3)
            return TestSetupFaultSSSE3KAT;
        if (backend == LEO2_BACKEND_AVX2)
            return TestSetupFaultAVX2KAT;
        if (backend == LEO2_BACKEND_AVX512)
            return TestSetupFaultAVX512KAT;
        if (backend == LEO2_BACKEND_GFNI)
            return TestSetupFaultGFNIKAT;
        throw std::runtime_error("unknown failure-backend argument");
    }
    throw std::runtime_error("unknown failure-stage argument");
}

bool runtime_can_select(leo2_backend backend, leo2_backend default_backend)
{
    (void)default_backend;
    return leopard::backend::TestBackendCanQualifyForHost(backend);
}

void check_table_accounting(
    leo2_backend backend,
    const leopard::backend::TestBackendState& state)
{
    if (backend == LEO2_BACKEND_SCALAR)
    {
        require(state.ff8_bytes == 65536U,
            "scalar GF8 table byte accounting changed");
        require(state.ff16_bytes == 8388608U,
            "scalar GF16 table byte accounting changed");
    }
    else if (backend == LEO2_BACKEND_GFNI)
    {
        // The GFNI member packs the four 8x8 affine blocks into 32 bytes per
        // GF16 logarithm (2 MiB) instead of the nibble shape's 128 bytes;
        // its GF8 table keeps the 32-byte duplicated-row shape.
        require(state.ff8_bytes == 8192U,
            "GFNI GF8 table byte accounting changed");
        require(state.ff16_bytes == 2097152U,
            "GFNI packed GF16 table byte accounting changed");
    }
    else
    {
        require(state.ff8_bytes == 8192U,
            "SIMD GF8 table byte accounting changed");
#if defined(LEO2_GFNI_VARIANT)
        // The documented global experiment retains AVX2/AVX-512 public
        // identities while replacing their GF16 nibble tables with the same
        // packed four-block affine representation as the production GFNI
        // member.  Scalar and SSSE3 remain separate translation units with
        // their ordinary table representations.
        if (backend == LEO2_BACKEND_AVX2 ||
            backend == LEO2_BACKEND_AVX512)
        {
            require(state.ff16_bytes == 2097152U,
                "global GFNI affine GF16 table byte accounting changed");
        }
        else
        {
            require(state.ff16_bytes == 8388608U,
                "non-GFNI SIMD GF16 table byte accounting changed");
        }
#else
        require(state.ff16_bytes == 8388608U,
            "SIMD GF16 table byte accounting changed");
#endif
    }
}

void run_failure_case(leo2_backend backend, const char* stage)
{
    using namespace leopard::backend;
    TestBackendState initial;
    if (!TestGetBackendState(backend, &initial))
    {
        std::printf("Backend failure case skipped: backend not compiled\n");
        return;
    }
    check_table_accounting(backend, initial);
    require(!initial.ff8_published && !initial.ff16_published &&
            !initial.qualified,
        "backend was unexpectedly initialized before fault injection");

    const leo2_backend default_backend = TestDefaultBackendForHost();
    require(default_backend != LEO2_BACKEND_AUTO,
        "host has no selectable process default");
    if (default_backend == LEO2_BACKEND_AVX512 &&
        backend == LEO2_BACKEND_AVX2 &&
        std::strcmp(stage, "kat") != 0)
    {
        // The AVX-512VL table deliberately reuses the AVX2 field tables.
        // In a forced-AVX512 process an AVX2 allocation fault is therefore a
        // startup dependency failure attributed to the selected AVX-512
        // backend, not a lazy AVX2 qualification failure.  AUTO/AVX2 builds
        // retain the dedicated AVX2 allocation-failure coverage.
        std::printf("Backend failure case skipped: AVX2 tables are an "
            "AVX-512 startup dependency\n");
        return;
    }
    if (!runtime_can_select(backend, default_backend))
    {
        std::printf("Backend failure case skipped: backend unavailable on host\n");
        return;
    }

    const bool kat = std::strcmp(stage, "kat") == 0;
    const leo2_result expected = kat ? LEO2_INTERNAL_ERROR
                                     : LEO2_OUT_OF_MEMORY;
    TestSetSetupFault(fault_for(backend, stage));

    static const unsigned kThreads = 24;
    std::atomic<unsigned> ready(0);
    std::atomic<bool> go(false);
    std::vector<leo2_result> results(kThreads, LEO2_SUCCESS);
    std::vector<leo2_context*> contexts(kThreads,
        reinterpret_cast<leo2_context*>(static_cast<uintptr_t>(1)));
    std::vector<std::thread> threads;
    for (unsigned i = 0; i < kThreads; ++i)
    {
        threads.push_back(std::thread([&, i]() {
            leo2_context_options options;
            std::memset(&options, 0, sizeof(options));
            options.struct_size = sizeof(options);
            options.backend = backend;
            options.thread_count = 1;
            ready.fetch_add(1, std::memory_order_release);
            while (!go.load(std::memory_order_acquire))
                std::this_thread::yield();
            results[i] = leo2_context_create(&options, &contexts[i]);
        }));
    }
    while (ready.load(std::memory_order_acquire) != kThreads)
        std::this_thread::yield();
    go.store(true, std::memory_order_release);
    for (size_t i = 0; i < threads.size(); ++i)
        threads[i].join();

    for (unsigned i = 0; i < kThreads; ++i)
    {
        require(results[i] == expected,
            "concurrent failure returned the wrong public result");
        require(contexts[i] == NULL,
            "failed context creation did not clear its output");
    }
    require(!TestSetupFaultPending(), "injected fault was not consumed");
    require(TestSetupFaultConsumptions() == 1,
        "one-shot setup fault was consumed more than once");

    leo2_context_options options;
    std::memset(&options, 0, sizeof(options));
    options.struct_size = sizeof(options);
    options.backend = backend;
    options.thread_count = 1;
    leo2_context* cached = reinterpret_cast<leo2_context*>(
        static_cast<uintptr_t>(1));
    require(leo2_context_create(&options, &cached) == expected,
        "cached failure changed after the injection was consumed");
    require(cached == NULL, "cached failure did not clear context output");
    require(TestSetupFaultConsumptions() == 1,
        "cached failure retried backend construction");

    TestBackendState failed;
    require(TestGetBackendState(backend, &failed),
        "failed backend state is unavailable");
    check_table_accounting(backend, failed);
    require(!failed.qualified,
        "failed backend was published as qualified");
    require(failed.failure == (kat ? QualificationSelfTestFailed
                                  : QualificationOutOfMemory),
        "cached internal qualification status is wrong");
    if (kat)
        require(failed.ff8_published && failed.ff16_published,
            "KAT failure did not retain one complete immutable table pair");
    else
        require(!failed.ff8_published && !failed.ff16_published,
            "allocation failure partially published backend tables");

    // A failed explicit lower-backend request must not mutate the already
    // qualified process default or AUTO selection.
    if (backend != default_backend)
    {
        leo2_context* automatic = NULL;
        leo2_context_options automatic_options;
        std::memset(&automatic_options, 0, sizeof(automatic_options));
        automatic_options.struct_size = sizeof(automatic_options);
        automatic_options.backend = LEO2_BACKEND_AUTO;
        automatic_options.thread_count = 1;
        require(leo2_context_create(&automatic_options, &automatic) ==
                LEO2_SUCCESS,
            "lower-backend failure changed AUTO availability");
        require(automatic != NULL &&
                leo2_context_backend(automatic) == default_backend &&
                SelectedBackend() == default_backend,
            "lower-backend failure changed the process default");
        leo2_context_destroy(automatic);
    }

    // Legacy initialization remains explicitly retryable.  Leopard2 must keep
    // the first context-setup outcome deterministic even if a direct legacy
    // retry subsequently succeeds and updates backend qualification globals.
    require(leo_init() == Leopard_Success,
        "direct legacy initialization retry failed");
    leo2_context* after_legacy_retry = reinterpret_cast<leo2_context*>(
        static_cast<uintptr_t>(1));
    require(leo2_context_create(&options, &after_legacy_retry) == expected,
        "legacy retry rewrote cached Leopard2 startup failure");
    require(after_legacy_retry == NULL,
        "cached failure after legacy retry did not clear context output");

    std::printf("Backend failure case passed: backend=%u stage=%s "
        "threads=%u expected=%d table_bytes=%llu\n",
        static_cast<unsigned>(backend), stage, kThreads,
        static_cast<int>(expected),
        static_cast<unsigned long long>(failed.ff8_bytes + failed.ff16_bytes));
}

void run_auto_avx512_fallback_case()
{
    using namespace leopard::backend;
    if (!IsCalibratedAutoAVX512EncodeHost() ||
        TestDefaultBackendForHost() != LEO2_BACKEND_AVX2 ||
        !TestBackendCanQualifyForHost(LEO2_BACKEND_AVX512))
    {
        std::printf("AUTO AVX-512 fallback skipped: host is not calibrated\n");
        return;
    }

    TestBackendState initial;
    require(TestGetBackendState(LEO2_BACKEND_AVX512, &initial),
        "AVX-512 backend state is unavailable");
    require(!initial.qualified,
        "AVX-512 unexpectedly qualified before fallback injection");
    TestSetSetupFault(TestSetupFaultAVX512KAT);

    leo2_context_options context_options;
    std::memset(&context_options, 0, sizeof(context_options));
    context_options.struct_size = sizeof(context_options);
    context_options.backend = LEO2_BACKEND_AUTO;
    context_options.thread_count = 1;
    leo2_context* context = NULL;
    require(leo2_context_create(&context_options, &context) == LEO2_SUCCESS &&
            context && leo2_context_backend(context) == LEO2_BACKEND_AVX2,
        "AUTO baseline context creation failed");

    leo2_codec* codec = NULL;
    require(leo2_codec_create(context, 8, 8,
            LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF16, NULL, &codec) ==
            LEO2_SUCCESS && codec,
        "optional AVX-512 KAT failure escaped codec setup");
    require(!TestSetupFaultPending() && TestSetupFaultConsumptions() == 1,
        "optional AVX-512 KAT fault consumption is wrong");
    leo2_backend selected = LEO2_BACKEND_AUTO;
    require(leo2_test_codec_transform_encode_backend(
            codec, 64, 8, 8, &selected) == LEO2_SUCCESS &&
            selected == LEO2_BACKEND_AVX2,
        "AUTO did not retain AVX2 after optional AVX-512 KAT failure");

    std::vector<std::vector<uint8_t> > original(
        8, std::vector<uint8_t>(64));
    std::vector<std::vector<uint8_t> > recovery(
        8, std::vector<uint8_t>(64));
    std::vector<const void*> original_ptrs(8);
    std::vector<void*> recovery_ptrs(8);
    for (unsigned shard = 0; shard < 8; ++shard)
    {
        for (unsigned i = 0; i < 64; ++i)
            original[shard][i] = static_cast<uint8_t>(
                shard * 37U + i * 19U);
        original_ptrs[shard] = original[shard].data();
        recovery_ptrs[shard] = recovery[shard].data();
    }
    size_t scratch_bytes = 0;
    require(leo2_encode_scratch_size(codec, 64, &scratch_bytes) ==
            LEO2_SUCCESS && scratch_bytes != 0,
        "fallback scratch query failed");
    std::vector<uint8_t> scratch_storage(
        scratch_bytes + leo2_scratch_alignment());
    const uintptr_t unaligned = reinterpret_cast<uintptr_t>(
        scratch_storage.data());
    const uintptr_t aligned = (unaligned + leo2_scratch_alignment() - 1U) &
        ~(static_cast<uintptr_t>(leo2_scratch_alignment()) - 1U);
    require(leo2_encode(codec, 64, original_ptrs.data(),
            recovery_ptrs.data(), reinterpret_cast<void*>(aligned),
            scratch_bytes) == LEO2_SUCCESS,
        "AVX2 fallback encode failed");

    leo2_codec* second = NULL;
    require(leo2_codec_create(context, 16, 16,
            LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF16, NULL, &second) ==
            LEO2_SUCCESS && second,
        "cached optional failure changed later codec setup");
    require(TestSetupFaultConsumptions() == 1,
        "cached optional failure retried AVX-512 qualification");
    leo2_codec_destroy(second);
    leo2_codec_destroy(codec);
    leo2_context_destroy(context);
    std::printf("AUTO AVX-512 KAT fallback passed\n");
}

void run_auto_gfni_encode_fallback_case(
    leopard::backend::TestSetupFault fault,
    const char* label)
{
    using namespace leopard::backend;
    if (!IsCalibratedAutoGF16GFNIEncodeHost() ||
        TestDefaultBackendForHost() != LEO2_BACKEND_AVX2 ||
        !TestBackendCanQualifyForHost(LEO2_BACKEND_GFNI))
    {
        std::printf("AUTO GF16 GFNI encode fallback skipped: host is not "
            "calibrated\n");
        return;
    }

    TestBackendState initial;
    require(TestGetBackendState(LEO2_BACKEND_GFNI, &initial),
        "GFNI backend state is unavailable");
    require(!initial.qualified,
        "GFNI unexpectedly qualified before fallback injection");
    require(leopard2_internal::SetAutoGF16GFNIEncodeEnabledForDiagnostics(
            true),
        "enable AUTO GF16 GFNI encode fallback probe");
    TestSetSetupFault(fault);

    leo2_context_options context_options;
    std::memset(&context_options, 0, sizeof(context_options));
    context_options.struct_size = sizeof(context_options);
    context_options.backend = LEO2_BACKEND_AUTO;
    context_options.thread_count = 1;
    leo2_context* context = NULL;
    require(leo2_context_create(&context_options, &context) == LEO2_SUCCESS &&
            context && leo2_context_backend(context) == LEO2_BACKEND_AVX2,
        "AUTO GF16 baseline context creation failed");

    leo2_codec* codec = NULL;
    require(leo2_codec_create(context, 1000, 200,
            LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF16, NULL, &codec) ==
            LEO2_SUCCESS && codec,
        "optional GFNI qualification failure escaped codec setup");
    require(!TestSetupFaultPending() && TestSetupFaultConsumptions() == 1,
        "optional GFNI qualification fault consumption is wrong");
    require(!leopard2_internal::AutoGF16GFNIEncodeAvailableForDiagnostics(
            codec),
        "failed GFNI table remained available to AUTO encode");
    leo2_backend selected = LEO2_BACKEND_AUTO;
    require(leo2_test_codec_transform_encode_backend(
            codec, 64U * 1024U, 200, 200, &selected) == LEO2_SUCCESS &&
            selected == LEO2_BACKEND_AVX2,
        "AUTO did not retain AVX2 after optional GFNI failure");

    static const size_t kBytes = 64U * 1024U;
    // Public encode permits source/source aliasing.  Reuse one source shard so
    // the exact eligible fallback call proves routing without allocating a
    // redundant 64-MiB input fixture.
    std::vector<uint8_t> original(kBytes);
    for (size_t i = 0; i < kBytes; ++i)
        original[i] = static_cast<uint8_t>(i * 19U);
    std::vector<uint8_t> recovery(200U * kBytes);
    std::vector<const void*> original_ptrs(1000, original.data());
    std::vector<void*> recovery_ptrs(200);
    for (unsigned shard = 0; shard < 200; ++shard)
        recovery_ptrs[shard] = recovery.data() + shard * kBytes;
    size_t scratch_bytes = 0;
    require(leo2_encode_scratch_size(codec, kBytes, &scratch_bytes) ==
            LEO2_SUCCESS && scratch_bytes != 0,
        "GFNI fallback scratch query failed");
    std::vector<uint8_t> scratch_storage(
        scratch_bytes + leo2_scratch_alignment());
    const uintptr_t unaligned = reinterpret_cast<uintptr_t>(
        scratch_storage.data());
    const uintptr_t aligned =
        (unaligned + leo2_scratch_alignment() - 1U) &
        ~(static_cast<uintptr_t>(leo2_scratch_alignment()) - 1U);
    require(leo2_encode(codec, kBytes, original_ptrs.data(),
            recovery_ptrs.data(), reinterpret_cast<void*>(aligned),
            scratch_bytes) == LEO2_SUCCESS,
        "AVX2 encode after optional GFNI failure failed");
    require(leopard2_internal::AutoGF16GFNIEncodeCallCountForDiagnostics() ==
            0U,
        "failed GFNI route counted an encode call");

    leo2_codec* second = NULL;
    require(leo2_codec_create(context, 1000, 200,
            LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF16, NULL, &second) ==
            LEO2_SUCCESS && second,
        "cached optional GFNI failure changed later codec setup");
    require(TestSetupFaultConsumptions() == 1,
        "cached optional GFNI failure retried qualification");
    require(!leopard2_internal::AutoGF16GFNIEncodeAvailableForDiagnostics(
            second),
        "cached failed GFNI table became available");
    require(leopard2_internal::FinishAutoGF16GFNIEncodeRouteProbeForDiagnostics(),
        "finish AUTO GF16 GFNI fallback probe");

    leo2_codec_destroy(second);
    leo2_codec_destroy(codec);
    leo2_context_destroy(context);
    std::printf("AUTO GF16 GFNI encode fallback passed: %s\n", label);
}

void run_auto_gfni_encode_disabled_inert_case()
{
    using namespace leopard::backend;
    if (!IsCalibratedAutoGF16GFNIEncodeHost() ||
        TestDefaultBackendForHost() != LEO2_BACKEND_AVX2 ||
        !TestBackendCanQualifyForHost(LEO2_BACKEND_GFNI))
    {
        std::printf("AUTO GF16 GFNI disabled-inert test skipped: host is "
            "not calibrated\n");
        return;
    }

    TestBackendState initial;
    require(TestGetBackendState(LEO2_BACKEND_GFNI, &initial),
        "GFNI backend state is unavailable");
    require(!initial.qualified && !initial.ff8_published &&
            !initial.ff16_published,
        "GFNI unexpectedly initialized before disabled-inert test");
    require(leopard2_internal::AutoGF16GFNIEncodeModeForDiagnostics() == 2U,
        "AUTO GF16 GFNI encode production default is not disabled");
    TestSetSetupFault(TestSetupFaultGFNIKAT);

    leo2_context_options context_options;
    std::memset(&context_options, 0, sizeof(context_options));
    context_options.struct_size = sizeof(context_options);
    context_options.backend = LEO2_BACKEND_AUTO;
    context_options.thread_count = 1;
    leo2_context* context = NULL;
    require(leo2_context_create(&context_options, &context) == LEO2_SUCCESS &&
            context && leo2_context_backend(context) == LEO2_BACKEND_AVX2,
        "disabled-inert AUTO context creation failed");

    leo2_codec* codec = NULL;
    require(leo2_codec_create(context, 1000, 200,
            LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF16, NULL, &codec) ==
            LEO2_SUCCESS && codec,
        "disabled-inert exact codec creation failed");
    require(TestSetupFaultPending() && TestSetupFaultConsumptions() == 0,
        "disabled selector initialized or qualified GFNI");
    require(!leopard2_internal::AutoGF16GFNIEncodeAvailableForDiagnostics(
                codec),
        "disabled selector retained a GFNI table");
    leo2_backend selected = LEO2_BACKEND_AUTO;
    require(leo2_test_codec_transform_encode_backend(
            codec, 64U * 1024U, 200, 200, &selected) == LEO2_SUCCESS &&
            selected == LEO2_BACKEND_AVX2,
        "disabled selector changed the exact encode backend");

    static const size_t kBytes = 64U * 1024U;
    std::vector<uint8_t> original(kBytes, 0x5a);
    std::vector<uint8_t> recovery(200U * kBytes);
    std::vector<const void*> original_ptrs(1000, original.data());
    std::vector<void*> recovery_ptrs(200);
    for (unsigned shard = 0; shard < 200; ++shard)
        recovery_ptrs[shard] = recovery.data() + shard * kBytes;
    size_t scratch_bytes = 0;
    require(leo2_encode_scratch_size(codec, kBytes, &scratch_bytes) ==
            LEO2_SUCCESS && scratch_bytes != 0,
        "disabled-inert scratch query failed");
    std::vector<uint8_t> scratch_storage(
        scratch_bytes + leo2_scratch_alignment());
    const uintptr_t unaligned = reinterpret_cast<uintptr_t>(
        scratch_storage.data());
    const uintptr_t aligned =
        (unaligned + leo2_scratch_alignment() - 1U) &
        ~(static_cast<uintptr_t>(leo2_scratch_alignment()) - 1U);
    require(leo2_encode(codec, kBytes, original_ptrs.data(),
            recovery_ptrs.data(), reinterpret_cast<void*>(aligned),
            scratch_bytes) == LEO2_SUCCESS,
        "disabled-inert AVX2 encode failed");
    require(TestSetupFaultPending() && TestSetupFaultConsumptions() == 0,
        "disabled exact encode touched GFNI qualification");
    TestBackendState final_state;
    require(TestGetBackendState(LEO2_BACKEND_GFNI, &final_state) &&
            !final_state.qualified && !final_state.ff8_published &&
            !final_state.ff16_published,
        "disabled exact encode allocated or published GFNI tables");
    require(leopard2_internal::AutoGF16GFNIEncodeCallCountForDiagnostics() ==
            0U,
        "disabled exact encode counted a GFNI route");
    require(leopard2_internal::AutoGF16GFNIEncodeModeForDiagnostics() == 2U,
        "disabled exact encode changed the production mode");

    TestSetSetupFault(TestSetupFaultNone);
    leo2_codec_destroy(codec);
    leo2_context_destroy(context);
    std::printf("AUTO GF16 GFNI encode disabled-inert test passed\n");
}

} // namespace

int main(int argc, char** argv)
{
    try
    {
        if (argc == 2 &&
            std::strcmp(argv[1], "auto-avx512-kat-fallback") == 0)
        {
            run_auto_avx512_fallback_case();
            return 0;
        }
        if (argc == 2 &&
            std::strcmp(argv[1], "auto-gfni-encode-disabled-inert") == 0)
        {
            run_auto_gfni_encode_disabled_inert_case();
            return 0;
        }
        if (argc == 2 &&
            std::strcmp(argv[1], "auto-gfni-encode-kat-fallback") == 0)
        {
            run_auto_gfni_encode_fallback_case(
                leopard::backend::TestSetupFaultGFNIKAT, "kat");
            return 0;
        }
        if (argc == 2 &&
            std::strcmp(argv[1],
                "auto-gfni-encode-ff8-allocation-fallback") == 0)
        {
            run_auto_gfni_encode_fallback_case(
                leopard::backend::TestSetupFaultGFNIFF8Allocation,
                "ff8-allocation");
            return 0;
        }
        if (argc == 2 &&
            std::strcmp(argv[1],
                "auto-gfni-encode-ff16-allocation-fallback") == 0)
        {
            run_auto_gfni_encode_fallback_case(
                leopard::backend::TestSetupFaultGFNIFF16Allocation,
                "ff16-allocation");
            return 0;
        }
        require(argc == 3,
            "usage: leopard2_backend_failures_test BACKEND STAGE");
        run_failure_case(parse_backend(argv[1]), argv[2]);
        return 0;
    }
    catch (const std::exception& error)
    {
        std::fprintf(stderr, "Leopard2 backend failure test failed: %s\n",
            error.what());
        return 1;
    }
}
