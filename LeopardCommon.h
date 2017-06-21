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

    Mid-term:
    + Add compile-time selectable XOR-only rowops instead of MULADD
    + Look into 12-bit fields as a performance optimization

    Long-term:
    + Evaluate the error locator polynomial based on fast polynomial interpolations in O(k log^2 k)
    + Look into getting EncodeL working so we can support larger recovery sets
    + Implement the decoder algorithm from {3} based on the Forney algorithm
*/

/*
    FFT Data Layout:

    We pack the data into memory in this order:

    [Recovery Data (Power of Two = M)] [Original Data] [Zero Padding out to 65536]

    For encoding, the placement is implied instead of actual memory layout.
    For decoding, the layout is explicitly used.
*/

/*
    Encoder algorithm:

    The encoder is described in {3}.  Operations are done O(K Log M),
    where K is the original data size, and M is up to twice the
    size of the recovery set.

    Roughly in brief:

        Recovery = FFT( IFFT(Data_0) xor IFFT(Data_1) xor ... )

    It walks the original data M chunks at a time performing the IFFT.
    Each IFFT intermediate result is XORed together into the first M chunks of
    the data layout.  Finally the FFT is performed.

    Encoder optimizations:
    * The first IFFT can be performed directly in the first M chunks.
    * The zero padding can be skipped while performing the final IFFT.
    Unrolling is used in the code to accomplish both these optimizations.
    * The final FFT can be truncated also if recovery set is not a power of 2.
    It is easy to truncate the FFT by ending the inner loop early.
    * The FFT operations can be unrolled two layers at a time so that instead
    of writing the result of the first layer out and reading it back in for
    the second layer, those interactions can happen in registers immediately.
*/

/*
    Decoder algorithm:

    The decoder is described in {1}.  Operations are done O(N Log N), where N is up
    to twice the size of the original data as described below.

    Roughly in brief:

        Original = -ErrLocator * FFT( Derivative( IFFT( ErrLocator * ReceivedData ) ) )


    Precalculations:
    ---------------

    At startup initialization, FFTInitialize() precalculates FWT(L) as
    described by equation (92) in {1}, where L = Log[i] for i = 0..Order,
    Order = 256 or 65536 for FF8/16.  This is stored in the LogWalsh vector.

    It also precalculates the FFT skew factors (s_i) as described by
    equation (28).  This is stored in the FFTSkew vector.

    For memory workspace N data chunks are needed, where N is a power of two
    at or above M + K.  K is the original data size and M is the next power
    of two above the recovery data size.  For example for K = 200 pieces of
    data and 10% redundancy, there are 20 redundant pieces, which rounds up
    to 32 = M.  M + K = 232 pieces, so N rounds up to 256.


    Online calculations:
    -------------------

    At runtime, the error locator polynomial is evaluated using the
    Fast Walsh-Hadamard transform as described in {1} equation (92).

    At runtime the data is explicit laid out in workspace memory like this:
    [Recovery Data (Power of Two = M)] [Original Data (K)] [Zero Padding out to N]

    Data that was lost is replaced with zeroes.
    Data that was received, including recovery data, is multiplied by the error
    locator polynomial as it is copied into the workspace.

    The IFFT is applied to the entire workspace of N chunks.
    Since the IFFT starts with pairs of inputs and doubles in width at each
    iteration, the IFFT is optimized by skipping zero padding at the end until
    it starts mixing with non-zero data.

    The formal derivative is applied to the entire workspace of N chunks.
    This is a massive XOR loop that runs 4 columns in parallel for speed.

    The FFT is applied to the entire workspace of N chunks.
    The FFT is optimized by only performing intermediate calculations required
    to recover lost data.  Since it starts wide and ends up working on adjacent
    pairs, at some point the intermediate results are not needed for data that
    will not be read by the application.  This optimization is implemented by
    the ErrorBitfield class.

    Finally, only recovered data is multiplied by the negative of the
    error locator polynomial as it is copied into the front of the
    workspace for the application to retrieve.
*/

/*
    Finite field arithmetic optimizations:

    For faster finite field multiplication, large tables are precomputed and
    applied during encoding/decoding on 64 bytes of data at a time using
    SSSE3 or AVX2 vector instructions and the ALTMAP approach from Jerasure.

    Addition in this finite field is XOR, and a vectorized memory XOR routine
    is also used.
*/

#include "leopard.h"

#include <stdint.h>
#ifdef _WIN32
#include <malloc.h>
#endif //_WIN32
#include <vector>
#include <atomic>
#include <memory>
#include <mutex>
#include <condition_variable>


//------------------------------------------------------------------------------
// Constants

// Enable 8-bit or 16-bit fields
#define LEO_HAS_FF8
#define LEO_HAS_FF16

// Enable using SIMD instructions
#define LEO_USE_SSSE3_OPT
#define LEO_USE_AVX2_OPT

// Avoid calculating final FFT values in decoder using bitfield
#define LEO_ERROR_BITFIELD_OPT

// Interleave butterfly operations between layer pairs in FFT
#define LEO_INTERLEAVE_BUTTERFLY4_OPT

// Optimize M=1 case
#define LEO_M1_OPT

// Unroll inner loops 4 times
#define LEO_USE_VECTOR4_OPT


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
// Windows Header

#ifdef _WIN32
    #define WIN32_LEAN_AND_MEAN

    #ifndef _WINSOCKAPI_
        #define DID_DEFINE_WINSOCKAPI
        #define _WINSOCKAPI_
    #endif
    #ifndef NOMINMAX
        #define NOMINMAX
    #endif
    #ifndef _WIN32_WINNT
        #define _WIN32_WINNT 0x0601 /* Windows 7+ */
    #endif

    #include <windows.h>
#endif

#ifdef DID_DEFINE_WINSOCKAPI
    #undef _WINSOCKAPI_
    #undef DID_DEFINE_WINSOCKAPI
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
    static_assert(sizeof(unsigned) == 4, "Assuming 32 bit unsigneds in LastNonzeroBit32");
    return 31 - (unsigned)__builtin_clz(x);
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

#ifdef LEO_M1_OPT

// x[] ^= y[] ^ z[]
void xor_mem_2to1(
    void * LEO_RESTRICT x,
    const void * LEO_RESTRICT y,
    const void * LEO_RESTRICT z,
    uint64_t bytes);

#endif // LEO_M1_OPT

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

// x[] ^= y[] (Multithreaded)
void VectorXOR_Threads(
    const uint64_t bytes,
    unsigned count,
    void** x,
    void** y);


//------------------------------------------------------------------------------
// XORSummer

class XORSummer
{
public:
    // Set the addition destination and byte count
    LEO_FORCE_INLINE void Initialize(void* dest)
    {
        DestBuffer = dest;
        Waiting = nullptr;
    }

    // Accumulate some source data
    LEO_FORCE_INLINE void Add(const void* src, const uint64_t bytes)
    {
#ifdef LEO_M1_OPT
        if (Waiting)
        {
            xor_mem_2to1(DestBuffer, src, Waiting, bytes);
            Waiting = nullptr;
        }
        else
            Waiting = src;
#else // LEO_M1_OPT
        xor_mem(DestBuffer, src, bytes);
#endif // LEO_M1_OPT
    }

    // Finalize in the destination buffer
    LEO_FORCE_INLINE void Finalize(const uint64_t bytes)
    {
#ifdef LEO_M1_OPT
        if (Waiting)
            xor_mem(DestBuffer, Waiting, bytes);
#endif // LEO_M1_OPT
    }

protected:
    void* DestBuffer;
    const void* Waiting;
};


//------------------------------------------------------------------------------
// SIMD-Safe Aligned Memory Allocations

static const unsigned kAlignmentBytes = LEO_ALIGN_BYTES;

static LEO_FORCE_INLINE uint8_t* SIMDSafeAllocate(size_t size)
{
    uint8_t* data = (uint8_t*)calloc(1, kAlignmentBytes + size);
    if (!data)
        return nullptr;
    unsigned offset = (unsigned)((uintptr_t)data % kAlignmentBytes);
    data += kAlignmentBytes - offset;
    data[-1] = (uint8_t)offset;
    return data;
}

static LEO_FORCE_INLINE void SIMDSafeFree(void* ptr)
{
    if (!ptr)
        return;
    uint8_t* data = (uint8_t*)ptr;
    unsigned offset = data[-1];
    if (offset >= kAlignmentBytes)
    {
        LEO_DEBUG_BREAK; // Should never happen
        return;
    }
    data -= kAlignmentBytes - offset;
    free(data);
}


} // namespace leopard
