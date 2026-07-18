/*
    Copyright (c) 2017 Christopher A. Taylor.  All rights reserved.
    See leopard.h for the full BSD license notice.
*/

#include "leopard.h"
#include "leopard2.h"

#include <atomic>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <thread>
#include <vector>

namespace {

static const unsigned kThreadCount = 64;

static void RecordFailure(
    std::atomic<unsigned>& failures,
    const char* operation,
    int result)
{
    std::fprintf(stderr, "%s failed: %d\n", operation, result);
    failures.fetch_add(1, std::memory_order_relaxed);
}

static void RunLegacyInitializer(
    unsigned index,
    std::atomic<unsigned>& ready,
    std::atomic<bool>& start,
    std::atomic<unsigned>& failures)
{
    ready.fetch_add(1, std::memory_order_release);
    while (!start.load(std::memory_order_acquire))
        std::this_thread::yield();

    // A version mismatch must never publish or poison initialization state.
    if ((index & 3u) == 0 &&
        leo_init_(LEO_VERSION + 1) != Leopard_InvalidInput)
        RecordFailure(failures, "invalid-version leo_init_", index);
    for (unsigned repeat = 0; repeat < 3; ++repeat)
    {
        const int result = leo_init();
        if (result != Leopard_Success)
            RecordFailure(failures, "legacy leo_init", result);
    }
}

static void RunContextInitializer(
    std::atomic<unsigned>& ready,
    std::atomic<bool>& start,
    std::atomic<unsigned>& failures)
{
    ready.fetch_add(1, std::memory_order_release);
    while (!start.load(std::memory_order_acquire))
        std::this_thread::yield();

    leo2_context_options options;
    std::memset(&options, 0, sizeof(options));
    options.struct_size = sizeof(options);
    options.thread_count = 1;
    for (unsigned repeat = 0; repeat < 2; ++repeat)
    {
        leo2_context* context = NULL;
        const leo2_result result = leo2_context_create(&options, &context);
        if (result != LEO2_SUCCESS || !context)
            RecordFailure(failures, "leo2_context_create", result);
        leo2_context_destroy(context);
    }
}

} // namespace

int main()
{
    // This executable intentionally reaches the library before any successful
    // initialization.  Invalid versions exercise retry/failure semantics.
    if (leo_init_(LEO_VERSION + 1) != Leopard_InvalidInput)
    {
        std::fprintf(stderr, "initial invalid-version result was not stable\n");
        return 1;
    }

    std::atomic<unsigned> ready(0);
    std::atomic<bool> start(false);
    std::atomic<unsigned> failures(0);
    std::vector<std::thread> threads;
    threads.reserve(kThreadCount);
    for (unsigned i = 0; i < kThreadCount; ++i)
    {
        if ((i & 1u) == 0)
            threads.push_back(std::thread(
                RunLegacyInitializer, i, std::ref(ready), std::ref(start),
                std::ref(failures)));
        else
            threads.push_back(std::thread(
                RunContextInitializer, std::ref(ready), std::ref(start),
                std::ref(failures)));
    }
    while (ready.load(std::memory_order_acquire) != kThreadCount)
        std::this_thread::yield();
    start.store(true, std::memory_order_release);
    for (size_t i = 0; i < threads.size(); ++i)
        threads[i].join();

    if (failures.load(std::memory_order_relaxed) != 0)
        return 2;
    if (leo_init() != Leopard_Success || leo_init() != Leopard_Success)
    {
        std::fprintf(stderr, "post-contention leo_init was not idempotent\n");
        return 3;
    }

    // Verify that the state published by the contending initializers is usable
    // by an ordinary legacy operation, not merely reported as successful.
    uint8_t original_storage[64];
    uint8_t recovery_storage[64];
    for (unsigned i = 0; i < sizeof(original_storage); ++i)
        original_storage[i] = static_cast<uint8_t>(i * 29u + 7u);
    const void* original[1] = { original_storage };
    void* work[1] = { recovery_storage };
    if (leo_encode(sizeof(original_storage), 1, 1, 1, original, work) !=
            Leopard_Success ||
        std::memcmp(original_storage, recovery_storage,
            sizeof(original_storage)) != 0)
    {
        std::fprintf(stderr, "published initialization state was unusable\n");
        return 4;
    }

    return 0;
}
