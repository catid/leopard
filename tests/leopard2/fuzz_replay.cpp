/* Deterministic smoke driver for the libFuzzer entry point. */

#include <stddef.h>
#include <stdint.h>
#include <errno.h>
#include <cstdlib>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef __has_feature
#define __has_feature(feature) 0
#endif

#if __has_feature(address_sanitizer) || defined(__SANITIZE_ADDRESS__)
#define LEO2_REPLAY_HAS_ASAN 1
#else
#define LEO2_REPLAY_HAS_ASAN 0
#endif

#if __has_feature(undefined_behavior_sanitizer)
#define LEO2_REPLAY_HAS_UBSAN 1
#else
#define LEO2_REPLAY_HAS_UBSAN 0
#endif

#ifndef LEO2_FUZZ_REPLAY_ROLE
#define LEO2_FUZZ_REPLAY_ROLE "unknown"
#endif

#if defined(__linux__) && (defined(__GNUC__) || defined(__clang__))
extern "C" int __asan_address_is_poisoned(const void*) __attribute__((weak));
// GCC treats the runtime handler's public spelling as a builtin and rejects
// taking its address.  Give the C++ identifier a private name while retaining
// an undefined weak reference to the exact ELF symbol used by compiler-rt.
extern "C" void leo2_ubsan_runtime_probe(void*, void*)
    __asm__("__ubsan_handle_type_mismatch_v1") __attribute__((weak));
extern "C" int __lsan_do_recoverable_leak_check() __attribute__((weak));
#define LEO2_REPLAY_ASAN_RUNTIME (__asan_address_is_poisoned != NULL)
#define LEO2_REPLAY_UBSAN_RUNTIME (leo2_ubsan_runtime_probe != NULL)
#define LEO2_REPLAY_LSAN_RUNTIME \
    (__lsan_do_recoverable_leak_check != NULL)
#else
#define LEO2_REPLAY_ASAN_RUNTIME 0
#define LEO2_REPLAY_UBSAN_RUNTIME 0
#define LEO2_REPLAY_LSAN_RUNTIME 0
#endif

extern "C" int LLVMFuzzerTestOneInput(const uint8_t* data, size_t size);
extern "C" int leo2_test_core_address_sanitizer_compiled();
extern "C" int leo2_test_core_undefined_sanitizer_compiled();
extern "C" size_t leo2_test_core_leak_sanitizer_canary();

#if defined(_MSC_VER)
#define LEO2_REPLAY_NOINLINE __declspec(noinline)
#elif defined(__GNUC__) || defined(__clang__)
#define LEO2_REPLAY_NOINLINE __attribute__((noinline))
#else
#define LEO2_REPLAY_NOINLINE
#endif

static LEO2_REPLAY_NOINLINE void scrub_leak_canary_roots()
{
    // A returned stack frame can retain a stale copy of an otherwise
    // unreachable allocation pointer.  Overwrite substantially more than the
    // tiny core canary frame before asking LSan for an explicit check.
    volatile uintptr_t scrub[2048];
    for (size_t i = 0; i < sizeof(scrub) / sizeof(scrub[0]); ++i)
        scrub[i] = 0;
}

static int recoverable_leak_check()
{
#if defined(__linux__) && (defined(__GNUC__) || defined(__clang__))
    if (!LEO2_REPLAY_LSAN_RUNTIME)
        return 0;
    return __lsan_do_recoverable_leak_check();
#else
    // Compiler-rt's LeakSanitizer interface is not available on this build.
    // The canary remains fail-closed through the caller's runtime check.
    return 0;
#endif
}

static uint64_t parse(const char* text)
{
    if (!text || !*text)
        return 0;
    errno = 0;
    char* end = NULL;
    const unsigned long long value = strtoull(text, &end, 0);
    return errno == 0 && end && *end == '\0' ? value : 0;
}

int main(int argc, char** argv)
{
    if (argc == 2 && strcmp(argv[1],
            "--leopard2-lsan-canary-v1") == 0)
    {
        const size_t leaked_bytes = leo2_test_core_leak_sanitizer_canary();
        scrub_leak_canary_roots();
        printf("leopard2_lsan_canary armed bytes=%llu\n",
            static_cast<unsigned long long>(leaked_bytes));
        fflush(stdout);
        if (leaked_bytes != 12345 || !LEO2_REPLAY_LSAN_RUNTIME)
            return 4;
        // The recoverable check emits compiler-rt's own diagnostic and returns
        // nonzero only when it found the now-unreachable linked-core object.
        // Returning a dedicated nonzero status makes a silent or disabled LSan
        // runtime fail closed rather than looking like a successful canary.
        if (recoverable_leak_check() <= 0)
            return 0;
        fflush(stderr);
        std::_Exit(86);
    }
    if (argc == 2 && strcmp(argv[1],
            "--leopard2-sanitizer-attestation-v1") == 0)
    {
        printf("{\"schema\":\"leopard2-sanitizer-attestation/v1\","
               "\"role\":\"%s\",\"driver\":\"deterministic-replay-v1\","
               "\"address_compile\":%s,\"address_runtime\":%s,"
               "\"undefined_compile\":%s,\"undefined_runtime\":%s,"
               "\"core_address_compile\":%s,"
               "\"core_undefined_compile\":%s}\n",
            LEO2_FUZZ_REPLAY_ROLE,
            LEO2_REPLAY_HAS_ASAN ? "true" : "false",
            LEO2_REPLAY_ASAN_RUNTIME ? "true" : "false",
            LEO2_REPLAY_HAS_UBSAN ? "true" : "false",
            LEO2_REPLAY_UBSAN_RUNTIME ? "true" : "false",
            leo2_test_core_address_sanitizer_compiled() ? "true" : "false",
            leo2_test_core_undefined_sanitizer_compiled() ? "true" : "false");
        return 0;
    }
    if (argc > 3)
        return 2;
    uint64_t state = argc > 1 ? parse(argv[1]) : UINT64_C(0x4c656f7061726432);
    const uint64_t iterations = argc > 2 ? parse(argv[2]) : 512;
    if (state == 0 || iterations == 0 || iterations > UINT64_C(10000000))
        return 2;
    uint8_t input[96];
    for (uint64_t iteration = 0; iteration < iterations; ++iteration)
    {
        for (size_t i = 0; i < sizeof(input); ++i)
        {
            state ^= state >> 12;
            state ^= state << 25;
            state ^= state >> 27;
            input[i] = static_cast<uint8_t>(
                state * UINT64_C(2685821657736338717));
        }
        if (LLVMFuzzerTestOneInput(input, sizeof(input)) != 0)
            return 1;
    }
    printf("leopard2_fuzz_replay seed=%llu iterations=%llu passed\n",
        static_cast<unsigned long long>(argc > 1 ? parse(argv[1]) :
            UINT64_C(0x4c656f7061726432)),
        static_cast<unsigned long long>(iterations));
    return 0;
}
