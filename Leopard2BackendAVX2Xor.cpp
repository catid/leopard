/*
    Copyright (c) 2017 Christopher A. Taylor.  All rights reserved.
    See Leopard2Backend.h for the full BSD license notice.

    Exact AVX2 XOR remainder arities live in this separate translation unit so
    they cannot change GCC's inlining or top-level ordering for the mature
    transform and reduction kernels in Leopard2BackendAVX2.cpp.
*/

#include "Leopard2Backend.h"

#include <cassert>
#include <cstring>
#include <immintrin.h>

namespace leopard { namespace backend {

#if !defined(NO_LEO_HAS_FF8)

#if defined(_MSC_VER)
#define LEO2_AVX2_XOR_FORCE_INLINE inline __forceinline
#elif defined(__GNUC__) || defined(__clang__)
#define LEO2_AVX2_XOR_FORCE_INLINE inline __attribute__((always_inline))
#else
#define LEO2_AVX2_XOR_FORCE_INLINE inline
#endif

#if defined(_MSC_VER)
#define LEO2_AVX2_FIXED_NOINLINE __declspec(noinline)
#elif defined(__GNUC__) && !defined(__clang__) && defined(__ELF__)
#define LEO2_AVX2_FIXED_NOINLINE \
    __attribute__((noinline, noclone, section(".text.leo2_avx2_r1_fixed"), \
        aligned(64)))
#elif (defined(__GNUC__) || defined(__clang__)) && defined(__ELF__)
#define LEO2_AVX2_FIXED_NOINLINE \
    __attribute__((noinline, section(".text.leo2_avx2_r1_fixed"), aligned(64)))
#elif defined(__GNUC__) || defined(__clang__)
#define LEO2_AVX2_FIXED_NOINLINE __attribute__((noinline))
#else
#define LEO2_AVX2_FIXED_NOINLINE
#endif

static LEO2_AVX2_XOR_FORCE_INLINE void AVX2AccumulateFixed64(
    const void* source_pointer,
    __m256i& accumulator0,
    __m256i& accumulator1)
{
    const uint8_t* const source =
        static_cast<const uint8_t*>(source_pointer);
    assert(source != NULL);
    accumulator0 = _mm256_xor_si256(accumulator0,
        _mm256_loadu_si256(reinterpret_cast<const __m256i*>(source)));
    accumulator1 = _mm256_xor_si256(accumulator1,
        _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(source + 32)));
}

static LEO2_AVX2_XOR_FORCE_INLINE void AVX2AccumulateFixed256(
    const void* source_pointer,
    __m256i& accumulator0,
    __m256i& accumulator1,
    __m256i& accumulator2,
    __m256i& accumulator3,
    __m256i& accumulator4,
    __m256i& accumulator5,
    __m256i& accumulator6,
    __m256i& accumulator7)
{
    const uint8_t* const source =
        static_cast<const uint8_t*>(source_pointer);
    assert(source != NULL);
    accumulator0 = _mm256_xor_si256(accumulator0,
        _mm256_loadu_si256(reinterpret_cast<const __m256i*>(source)));
    accumulator1 = _mm256_xor_si256(accumulator1,
        _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(source + 32)));
    accumulator2 = _mm256_xor_si256(accumulator2,
        _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(source + 64)));
    accumulator3 = _mm256_xor_si256(accumulator3,
        _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(source + 96)));
    accumulator4 = _mm256_xor_si256(accumulator4,
        _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(source + 128)));
    accumulator5 = _mm256_xor_si256(accumulator5,
        _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(source + 160)));
    accumulator6 = _mm256_xor_si256(accumulator6,
        _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(source + 192)));
    accumulator7 = _mm256_xor_si256(accumulator7,
        _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(source + 224)));
}

LEO2_AVX2_FIXED_NOINLINE void AVX2XorMemorySourcesFixed64(
    void* destination,
    const void* initial_source,
    const void* const* sources,
    uint32_t source_count,
    uint32_t excluded_source)
{
    assert(destination != NULL && initial_source != NULL && sources != NULL);
    assert(excluded_source <= source_count);
    assert(excluded_source == source_count ||
        sources[excluded_source] == NULL);
    const uint8_t* const initial =
        static_cast<const uint8_t*>(initial_source);
    __m256i accumulator0 = _mm256_loadu_si256(
        reinterpret_cast<const __m256i*>(initial));
    __m256i accumulator1 = _mm256_loadu_si256(
        reinterpret_cast<const __m256i*>(initial + 32));

    uint32_t source_index = 0;
    for (; source_index < excluded_source; ++source_index)
        AVX2AccumulateFixed64(
            sources[source_index], accumulator0, accumulator1);
    if (excluded_source < source_count)
        ++source_index;
    for (; source_index < source_count; ++source_index)
        AVX2AccumulateFixed64(
            sources[source_index], accumulator0, accumulator1);

    uint8_t* const output = static_cast<uint8_t*>(destination);
    _mm256_storeu_si256(
        reinterpret_cast<__m256i*>(output), accumulator0);
    _mm256_storeu_si256(
        reinterpret_cast<__m256i*>(output + 32), accumulator1);
}

LEO2_AVX2_FIXED_NOINLINE void AVX2XorMemorySourcesFixed256(
    void* destination,
    const void* initial_source,
    const void* const* sources,
    uint32_t source_count,
    uint32_t excluded_source)
{
    assert(destination != NULL && initial_source != NULL && sources != NULL);
    assert(excluded_source <= source_count);
    assert(excluded_source == source_count ||
        sources[excluded_source] == NULL);
    const uint8_t* const initial =
        static_cast<const uint8_t*>(initial_source);
    __m256i accumulator0 = _mm256_loadu_si256(
        reinterpret_cast<const __m256i*>(initial));
    __m256i accumulator1 = _mm256_loadu_si256(
        reinterpret_cast<const __m256i*>(initial + 32));
    __m256i accumulator2 = _mm256_loadu_si256(
        reinterpret_cast<const __m256i*>(initial + 64));
    __m256i accumulator3 = _mm256_loadu_si256(
        reinterpret_cast<const __m256i*>(initial + 96));
    __m256i accumulator4 = _mm256_loadu_si256(
        reinterpret_cast<const __m256i*>(initial + 128));
    __m256i accumulator5 = _mm256_loadu_si256(
        reinterpret_cast<const __m256i*>(initial + 160));
    __m256i accumulator6 = _mm256_loadu_si256(
        reinterpret_cast<const __m256i*>(initial + 192));
    __m256i accumulator7 = _mm256_loadu_si256(
        reinterpret_cast<const __m256i*>(initial + 224));

    uint32_t source_index = 0;
    for (; source_index < excluded_source; ++source_index)
        AVX2AccumulateFixed256(sources[source_index],
            accumulator0, accumulator1, accumulator2, accumulator3,
            accumulator4, accumulator5, accumulator6, accumulator7);
    if (excluded_source < source_count)
        ++source_index;
    for (; source_index < source_count; ++source_index)
        AVX2AccumulateFixed256(sources[source_index],
            accumulator0, accumulator1, accumulator2, accumulator3,
            accumulator4, accumulator5, accumulator6, accumulator7);

    uint8_t* const output = static_cast<uint8_t*>(destination);
    _mm256_storeu_si256(
        reinterpret_cast<__m256i*>(output), accumulator0);
    _mm256_storeu_si256(
        reinterpret_cast<__m256i*>(output + 32), accumulator1);
    _mm256_storeu_si256(
        reinterpret_cast<__m256i*>(output + 64), accumulator2);
    _mm256_storeu_si256(
        reinterpret_cast<__m256i*>(output + 96), accumulator3);
    _mm256_storeu_si256(
        reinterpret_cast<__m256i*>(output + 128), accumulator4);
    _mm256_storeu_si256(
        reinterpret_cast<__m256i*>(output + 160), accumulator5);
    _mm256_storeu_si256(
        reinterpret_cast<__m256i*>(output + 192), accumulator6);
    _mm256_storeu_si256(
        reinterpret_cast<__m256i*>(output + 224), accumulator7);
}

#undef LEO2_AVX2_FIXED_NOINLINE

template<unsigned SourceCount>
static LEO2_AVX2_XOR_FORCE_INLINE void AVX2XorMemoryFusedFinalGroup(
    void* destination,
    const void* initial_source,
    const void* const* sources,
    uint64_t byte_count)
{
    uint8_t* output = static_cast<uint8_t*>(destination);
    const uint8_t* initial = static_cast<const uint8_t*>(initial_source);
    const uint8_t* inputs[SourceCount];
    for (unsigned lane = 0; lane < SourceCount; ++lane)
        inputs[lane] = static_cast<const uint8_t*>(sources[lane]);

    uint64_t offset = 0;
    while (byte_count - offset >= 32)
    {
        __m256i result = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(initial + offset));
        for (unsigned lane = 0; lane < SourceCount; ++lane)
            result = _mm256_xor_si256(result, _mm256_loadu_si256(
                reinterpret_cast<const __m256i*>(inputs[lane] + offset)));
        _mm256_storeu_si256(
            reinterpret_cast<__m256i*>(output + offset), result);
        offset += 32;
    }

    /* The public API permits arbitrary positive byte lengths.  Volatile keeps
       GCC and Clang from cloning overlap-checking vector loops for this cold
       suffix; the API requires destination/source ranges to be disjoint. */
    volatile uint8_t* const tail_output = output;
    const volatile uint8_t* const tail_initial = initial;
    while (offset < byte_count)
    {
        uint8_t result = tail_initial[offset];
        for (unsigned lane = 0; lane < SourceCount; ++lane)
            result ^= static_cast<const volatile uint8_t*>(
                inputs[lane])[offset];
        tail_output[offset++] = result;
    }
}

void AVX2XorMemorySourcesFusedFinal(
    void* destination,
    const void* initial_source,
    const void* const* sources,
    uint32_t source_count,
    uint64_t byte_count)
{
    if (byte_count == 0)
        return;

    const void* waiting[8];
    unsigned waiting_count = 0;
    const void* accumulator = initial_source;
    bool wrote_destination = false;
    for (uint32_t i = 0; i < source_count; ++i)
    {
        if (!sources[i])
            continue;
        waiting[waiting_count++] = sources[i];
        if (waiting_count == 8)
        {
            AVX2XorMemoryFusedFinalGroup<8>(
                destination, accumulator, waiting, byte_count);
            accumulator = destination;
            wrote_destination = true;
            waiting_count = 0;
        }
    }

    switch (waiting_count)
    {
    case 7:
        AVX2XorMemoryFusedFinalGroup<7>(
            destination, accumulator, waiting, byte_count);
        return;
    case 6:
        AVX2XorMemoryFusedFinalGroup<6>(
            destination, accumulator, waiting, byte_count);
        return;
    case 5:
        AVX2XorMemoryFusedFinalGroup<5>(
            destination, accumulator, waiting, byte_count);
        return;
    case 4:
        AVX2XorMemoryFusedFinalGroup<4>(
            destination, accumulator, waiting, byte_count);
        return;
    case 3:
        AVX2XorMemoryFusedFinalGroup<3>(
            destination, accumulator, waiting, byte_count);
        return;
    case 2:
        AVX2XorMemoryFusedFinalGroup<2>(
            destination, accumulator, waiting, byte_count);
        return;
    case 1:
        AVX2XorMemoryFusedFinalGroup<1>(
            destination, accumulator, waiting, byte_count);
        return;
    default:
        break;
    }
    if (!wrote_destination)
        std::memcpy(destination, initial_source, static_cast<size_t>(byte_count));
}

void AVX2XorMemorySourcesGroup4(
    void* destination,
    const void* initial_source,
    const void* const* sources,
    uint32_t source_count,
    uint64_t byte_count)
{
    if (byte_count == 0)
        return;

    const void* waiting[4];
    unsigned waiting_count = 0;
    const void* accumulator = initial_source;
    bool wrote_destination = false;
    for (uint32_t i = 0; i < source_count; ++i)
    {
        if (!sources[i])
            continue;
        waiting[waiting_count++] = sources[i];
        if (waiting_count == 4)
        {
            AVX2XorMemoryFusedFinalGroup<4>(
                destination, accumulator, waiting, byte_count);
            accumulator = destination;
            wrote_destination = true;
            waiting_count = 0;
        }
    }

    switch (waiting_count)
    {
    case 3:
        AVX2XorMemoryFusedFinalGroup<3>(
            destination, accumulator, waiting, byte_count);
        return;
    case 2:
        AVX2XorMemoryFusedFinalGroup<2>(
            destination, accumulator, waiting, byte_count);
        return;
    case 1:
        AVX2XorMemoryFusedFinalGroup<1>(
            destination, accumulator, waiting, byte_count);
        return;
    default:
        break;
    }
    if (!wrote_destination)
        std::memcpy(destination, initial_source, static_cast<size_t>(byte_count));
}

#undef LEO2_AVX2_XOR_FORCE_INLINE

#endif // !NO_LEO_HAS_FF8

}} // namespace leopard::backend
