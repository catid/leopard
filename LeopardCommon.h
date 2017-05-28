/*
    Copyright (c) 2017 Christopher A. Taylor.  All rights reserved.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are met:

    * Redistributions of source code must retain the above copyright notice,
      this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright notice,
      this list of conditions and the following disclaimer in the documentation
      and/or other materials provided with the distribution.
    * Neither the name of Leopard-RS nor the names of its contributors may be
      used to endorse or promote products derived from this software without
      specific prior written permission.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
    AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
    IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
    ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
    LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
    CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
    SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
    INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
    CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
    ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
    POSSIBILITY OF SUCH DAMAGE.
*/

#pragma once

/*
    TODO:
    + Fixes for all different input sizes
    + New 16-bit Muladd inner loops
        + Benchmarks for large data!
    + Add multi-threading to split up long parallelizable calculations
        + Write detailed comments for all the routines
        + Final benchmarks!
    + Release version 1
        + Finish up documentation

    TBD:
    + Look into getting EncodeL working so we can support smaller data (Ask Lin)
    + Look into using FFT_m instead of FFT_n for decoder
*/

#include "leopard.h"

#include <stdint.h>


//------------------------------------------------------------------------------
// Constants

// Unroll inner loops 4 times
#define LEO_USE_VECTOR4_OPT

// Define this to enable the optimized version of FWHT()
#define LEO_FWHT_OPT

// Avoid scheduling reduced FFT operations that are unneeded
#define LEO_SCHEDULE_OPT


//------------------------------------------------------------------------------
// Debug

// Some bugs only repro in release mode, so this can be helpful
//#define LEO_DEBUG_IN_RELEASE

#if defined(_DEBUG) || defined(DEBUG) || defined(LEO_DEBUG_IN_RELEASE)
    #define LEO_DEBUG
    #ifdef _WIN32
        #define LEO_DEBUG_BREAK __debugbreak()
    #else
        #define LEO_DEBUG_BREAK __builtin_trap()
    #endif
    #define LEO_DEBUG_ASSERT(cond) { if (!(cond)) { LEO_DEBUG_BREAK; } }
#else
    #define LEO_DEBUG_BREAK ;
    #define LEO_DEBUG_ASSERT(cond) ;
#endif


//------------------------------------------------------------------------------
// Platform/Architecture

#ifdef _MSC_VER
    #include <intrin.h>
#endif

#if defined(ANDROID) || defined(IOS)
    #define LEO_TARGET_MOBILE
#endif // ANDROID

#if defined(__AVX2__) || (defined (_MSC_VER) && _MSC_VER >= 1900)
    #define LEO_TRY_AVX2 /* 256-bit */
    #include <immintrin.h>
    #define LEO_ALIGN_BYTES 32
#else // __AVX2__
    #define LEO_ALIGN_BYTES 16
#endif // __AVX2__

#if !defined(LEO_TARGET_MOBILE)
    // Note: MSVC currently only supports SSSE3 but not AVX2
    #include <tmmintrin.h> // SSSE3: _mm_shuffle_epi8
    #include <emmintrin.h> // SSE2
#endif // LEO_TARGET_MOBILE

#if defined(HAVE_ARM_NEON_H)
    #include <arm_neon.h>
#endif // HAVE_ARM_NEON_H

#if defined(LEO_TARGET_MOBILE)

    #define LEO_ALIGNED_ACCESSES /* Inputs must be aligned to LEO_ALIGN_BYTES */

# if defined(HAVE_ARM_NEON_H)
    // Compiler-specific 128-bit SIMD register keyword
    #define LEO_M128 uint8x16_t
    #define LEO_TRY_NEON
#else
    #define LEO_M128 uint64_t
# endif

#else // LEO_TARGET_MOBILE

    // Compiler-specific 128-bit SIMD register keyword
    #define LEO_M128 __m128i

#endif // LEO_TARGET_MOBILE

#ifdef LEO_TRY_AVX2
    // Compiler-specific 256-bit SIMD register keyword
    #define LEO_M256 __m256i
#endif

// Compiler-specific C++11 restrict keyword
#define LEO_RESTRICT __restrict

// Compiler-specific force inline keyword
#ifdef _MSC_VER
    #define LEO_FORCE_INLINE inline __forceinline
#else
    #define LEO_FORCE_INLINE inline __attribute__((always_inline))
#endif

// Compiler-specific alignment keyword
// Note: Alignment only matters for ARM NEON where it should be 16
#ifdef _MSC_VER
    #define LEO_ALIGNED __declspec(align(LEO_ALIGN_BYTES))
#else // _MSC_VER
    #define LEO_ALIGNED __attribute__((aligned(LEO_ALIGN_BYTES)))
#endif // _MSC_VER


namespace leopard {


//------------------------------------------------------------------------------
// Runtime CPU Architecture Check

// Initialize CPU architecture flags
void InitializeCPUArch();


#if defined(LEO_TRY_NEON)
# if defined(IOS) && defined(__ARM_NEON__)
    // Does device support NEON?
    static const bool CpuHasNeon = true;
    static const bool CpuHasNeon64 = true;
# else
    // Does device support NEON?
    // Remember to add LOCAL_STATIC_LIBRARIES := cpufeatures
    extern bool CpuHasNeon; // V6 / V7
    extern bool CpuHasNeon64; // 64-bit
# endif
#endif

#if !defined(LEO_TARGET_MOBILE)
# if defined(LEO_TRY_AVX2)
    // Does CPU support AVX2?
    extern bool CpuHasAVX2;
# endif
    // Does CPU support SSSE3?
    extern bool CpuHasSSSE3;
#endif // LEO_TARGET_MOBILE


//------------------------------------------------------------------------------
// Portable Intrinsics

// Returns highest bit index 0..31 where the first non-zero bit is found
// Precondition: x != 0
LEO_FORCE_INLINE unsigned LastNonzeroBit32(unsigned x)
{
#ifdef _MSC_VER
    unsigned long index;
    // Note: Ignoring result because x != 0
    _BitScanReverse(&index, (uint32_t)x);
    return (unsigned)index;
#else
    // Note: Ignoring return value of 0 because x != 0
    return 31 - (unsigned)__builtin_clzl(x);
#endif
}

// Returns next power of two at or above given value
LEO_FORCE_INLINE unsigned NextPow2(unsigned n)
{
    return 2UL << LastNonzeroBit32(n - 1);
}


//------------------------------------------------------------------------------
// XOR Memory
//
// This works for both 8-bit and 16-bit finite fields

// x[] ^= y[]
void xor_mem(
    void * LEO_RESTRICT x, const void * LEO_RESTRICT y,
    uint64_t bytes);

#ifdef LEO_USE_VECTOR4_OPT

// For i = {0, 1, 2, 3}: x_i[] ^= x_i[]
void xor_mem4(
    void * LEO_RESTRICT x_0, const void * LEO_RESTRICT y_0,
    void * LEO_RESTRICT x_1, const void * LEO_RESTRICT y_1,
    void * LEO_RESTRICT x_2, const void * LEO_RESTRICT y_2,
    void * LEO_RESTRICT x_3, const void * LEO_RESTRICT y_3,
    uint64_t bytes);

#endif // LEO_USE_VECTOR4_OPT

// x[] ^= y[]
void VectorXOR(
    const uint64_t bytes,
    unsigned count,
    void** x,
    void** y);


} // namespace leopard
