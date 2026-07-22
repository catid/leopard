/*
    Copyright (c) 2017 Christopher A. Taylor.  All rights reserved.
    See Leopard2Backend.h for the full BSD license notice.

    Exact AVX2 XOR remainder arities live in this separate translation unit so
    they cannot change GCC's inlining or top-level ordering for the mature
    transform and reduction kernels in Leopard2BackendAVX2.cpp.
*/

#include "Leopard2Backend.h"

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

#undef LEO2_AVX2_XOR_FORCE_INLINE

#endif // !NO_LEO_HAS_FF8

}} // namespace leopard::backend
