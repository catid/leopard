/*
    Copyright (c) 2017 Christopher A. Taylor.  All rights reserved.
    See Leopard2Backend.h for the full BSD license notice.
*/

#include "Leopard2Backend.h"
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
    throw std::runtime_error("unknown backend argument");
}

leopard::backend::TestSetupFault fault_for(
    leo2_backend backend,
    const char* stage)
{
    using namespace leopard::backend;
    if (std::strcmp(stage, "ff8-allocation") == 0)
    {
        if (backend == LEO2_BACKEND_SCALAR)
            return TestSetupFaultScalarFF8Allocation;
        if (backend == LEO2_BACKEND_SSSE3)
            return TestSetupFaultSSSE3FF8Allocation;
        if (backend == LEO2_BACKEND_AVX2)
            return TestSetupFaultAVX2FF8Allocation;
        return TestSetupFaultAVX512FF8Allocation;
    }
    if (std::strcmp(stage, "ff16-allocation") == 0)
    {
        if (backend == LEO2_BACKEND_SCALAR)
            return TestSetupFaultScalarFF16Allocation;
        if (backend == LEO2_BACKEND_SSSE3)
            return TestSetupFaultSSSE3FF16Allocation;
        if (backend == LEO2_BACKEND_AVX2)
            return TestSetupFaultAVX2FF16Allocation;
        return TestSetupFaultAVX512FF16Allocation;
    }
    if (std::strcmp(stage, "kat") == 0)
    {
        if (backend == LEO2_BACKEND_SCALAR)
            return TestSetupFaultScalarKAT;
        if (backend == LEO2_BACKEND_SSSE3)
            return TestSetupFaultSSSE3KAT;
        if (backend == LEO2_BACKEND_AVX2)
            return TestSetupFaultAVX2KAT;
        return TestSetupFaultAVX512KAT;
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
    else
    {
        require(state.ff8_bytes == 8192U,
            "SIMD GF8 table byte accounting changed");
        require(state.ff16_bytes == 8388608U,
            "SIMD GF16 table byte accounting changed");
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

} // namespace

int main(int argc, char** argv)
{
    try
    {
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
