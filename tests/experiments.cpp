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

#include <string.h>
#include <time.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>


//#define LEO_SHORT_FIELD
//#define LEO_EXPERIMENT_EXTRA_XOR
//#define LEO_EXPERIMENT_EXTRA_MULS
#define LEO_EXPERIMENT_CANTOR_BASIS

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

// Compiler-specific C++11 restrict keyword
#define LEO_RESTRICT __restrict

// Compiler-specific force inline keyword
#ifdef _MSC_VER
    #define LEO_FORCE_INLINE inline __forceinline
#else
    #define LEO_FORCE_INLINE inline __attribute__((always_inline))
#endif


//------------------------------------------------------------------------------
// Field

#ifdef LEO_SHORT_FIELD
typedef uint8_t ffe_t;
static const unsigned kGFBits = 8;
static const unsigned kGFPolynomial = 0x11D;
ffe_t kGFBasis[kGFBits] = {
#ifdef LEO_EXPERIMENT_CANTOR_BASIS
    1, 214, 152, 146, 86, 200, 88, 230 // Cantor basis
#else
    1, 2, 4, 8, 16, 32, 64, 128 // Monomial basis
#endif
};
#else
typedef uint16_t ffe_t;
static const unsigned kGFBits = 16;
static const unsigned kGFPolynomial = 0x1002D;
ffe_t kGFBasis[kGFBits] = {
#ifdef LEO_EXPERIMENT_CANTOR_BASIS
    0x0001, 0xACCA, 0x3C0E, 0x163E, // Cantor basis
    0xC582, 0xED2E, 0x914C, 0x4012,
    0x6C98, 0x10D8, 0x6A72, 0xB900,
    0xFDB8, 0xFB34, 0xFF38, 0x991E
#else
    1, 2, 4, 8, // Monomial basis
    16, 32, 64, 128,
    256, 512, 1024, 2048,
    4096, 8192, 16384, 32768
#endif
};
#endif

static const unsigned kFieldSize = (unsigned)1 << kGFBits; //Field size
static const unsigned kModulus = kFieldSize - 1;

static ffe_t GFLog[kFieldSize];
static ffe_t GFExp[kFieldSize];

// Initialize GFLog[], GFExp[]
static void InitField()
{
    unsigned state = 1;
    for (unsigned i = 0; i < kModulus; ++i)
    {
        GFExp[state] = static_cast<ffe_t>(i);
        state <<= 1;
        if (state >= kFieldSize)
            state ^= kGFPolynomial;
    }
    GFExp[0] = kModulus;

    // Conversion to chosen basis:

    GFLog[0] = 0;
    for (unsigned i = 0; i < kGFBits; ++i)
    {
        const ffe_t basis = kGFBasis[i];
        const unsigned width = (unsigned)(1UL << i);

        for (unsigned j = 0; j < width; ++j)
            GFLog[j + width] = GFLog[j] ^ basis;
    }

    for (unsigned i = 0; i < kFieldSize; ++i)
        GFLog[i] = GFExp[GFLog[i]];

    for (unsigned i = 0; i < kFieldSize; ++i)
        GFExp[GFLog[i]] = i;

    GFExp[kModulus] = GFExp[0];
}


//------------------------------------------------------------------------------
// Mod Q Field Operations
//
// Q is the maximum symbol value, e.g. 255 or 65535.

// z = x + y (mod Q)
static inline ffe_t AddModQ(ffe_t a, ffe_t b)
{
    const unsigned sum = (unsigned)a + b;

    // Partial reduction step, allowing for Q to be returned
    return static_cast<ffe_t>(sum + (sum >> kGFBits));
}

// z = x - y (mod Q)
static inline ffe_t SubModQ(ffe_t a, ffe_t b)
{
    const unsigned dif = (unsigned)a - b;

    // Partial reduction step, allowing for Q to be returned
    return static_cast<ffe_t>(dif + (dif >> kGFBits));
}

// return a*GFExp[b] over GF(2^r)
static ffe_t mulE(ffe_t a, ffe_t b)
{
    if (a == 0)
        return 0;
//    if (b == 0)
//        return a;

    const ffe_t sum = static_cast<ffe_t>(AddModQ(GFLog[a], b));
    return GFExp[sum];
}


//------------------------------------------------------------------------------
// Fast Walsh-Hadamard Transform (FWHT) Mod Q
//
// Q is the maximum symbol value, e.g. 255 or 65535.

// Define this to enable the optimized version of FWHT()
#define LEO_FWHT_OPTIMIZED

typedef ffe_t fwht_t;

// {a, b} = {a + b, a - b} (Mod Q)
static LEO_FORCE_INLINE void FWHT_2(fwht_t& LEO_RESTRICT a, fwht_t& LEO_RESTRICT b)
{
    const fwht_t sum = AddModQ(a, b);
    const fwht_t dif = SubModQ(a, b);
    a = sum;
    b = dif;
}

// Reference implementation
static void FWHT(fwht_t* data, const unsigned bits)
{
    const unsigned size = (unsigned)(1UL << bits);
    for (unsigned width = 1; width < size; width <<= 1)
        for (unsigned i = 0; i < size; i += (width << 1))
            for (unsigned j = i; j < (width + i); ++j)
                FWHT_2(data[j], data[j + width]);
}


//------------------------------------------------------------------------------
// Formal Derivative

// Formal derivative of polynomial in the new basis
static void formal_derivative(ffe_t* cos, const unsigned size)
{
    /*
        Left to right xoring data ahead into data behind.

        If the data ends in all zeroes, this can simply stop.
    */
    for (unsigned i = 1; i < size; ++i)
    {
        const unsigned leng = ((i ^ (i - 1)) + 1) >> 1;

        // If a large number of values are being XORed:
        for (unsigned j = i - leng; j < i; ++j)
            cos[j] ^= cos[j + leng];
    }

    // Doesn't seem to be needed
#ifdef LEO_EXPERIMENT_EXTRA_XOR
    /*
        Same here - Zeroes on the right are preserved
    */
    for (unsigned i = size; i < kFieldSize; i <<= 1)
    {
        for (unsigned j = 0; j < size; ++j)
            cos[j] ^= cos[j + i];
    }
#endif
}


//------------------------------------------------------------------------------
// Fast Fourier Transform

static ffe_t skewVec[kModulus]; // twisted factors used in FFT

static LEO_FORCE_INLINE void ifft_butterfly(ffe_t& a, ffe_t& b, ffe_t skew)
{
    b ^= a;
    if (skew != kModulus)
        a ^= mulE(b, skew);
}

// IFFT in the proposed basis
static void IFLT(ffe_t* data, const unsigned size, const unsigned index)
{
    for (unsigned width = 1; width < size; width <<= 1)
    {
        for (unsigned j = width; j < size; j += (width << 1))
        {
            const ffe_t skew = skewVec[j + index - 1];

            for (unsigned i = j - width; i < j; ++i)
                ifft_butterfly(data[i], data[i + width], skew);
        }
    }
}

static LEO_FORCE_INLINE void fft_butterfly(ffe_t& a, ffe_t& b, ffe_t skew)
{
    if (skew != kModulus)
        a ^= mulE(b, skew);
    b ^= a;
}

// FFT in the proposed basis
static void FLT(ffe_t* data, const unsigned size, const unsigned skewIndex, const unsigned output_elements)
{
    for (unsigned width = (size >> 1); width > 0; width >>= 1)
    {
        const ffe_t* skewLUT = skewVec + width + skewIndex - 1;

        for (unsigned j = 0; j < output_elements; j += (width << 1))
        {
            const ffe_t skew = skewLUT[j];

            for (unsigned i = j; i < j + width; ++i)
                fft_butterfly(data[i], data[i + width], skew);
        }
    }
}


//------------------------------------------------------------------------------
// FFT Initialization

#ifdef LEO_EXPERIMENT_EXTRA_MULS
static ffe_t B[kFieldSize >> 1];     // factors used in formal derivative
#endif
static fwht_t log_walsh[kFieldSize];  // factors used in the evaluation of the error locator polynomial

// Initialize skewVec[], B[], log_walsh[]
static void InitFieldOperations()
{
    ffe_t temp[kGFBits - 1];

    for (unsigned i = 1; i < kGFBits; ++i)
        temp[i - 1] = (ffe_t)((unsigned)1 << i);

    for (unsigned m = 0; m < (kGFBits - 1); ++m)
    {
        const unsigned step = (unsigned)1 << (m + 1);

        skewVec[((unsigned)1 << m) - 1] = 0;

        for (unsigned i = m; i < (kGFBits - 1); ++i)
        {
            const unsigned s = ((unsigned)1 << (i + 1));

            for (unsigned j = ((unsigned)1 << m) - 1; j < s; j += step)
                skewVec[j + s] = skewVec[j] ^ temp[i];
        }

        temp[m] = kModulus - GFLog[mulE(temp[m], GFLog[temp[m] ^ 1])];

        for (unsigned i = m + 1; i < (kGFBits - 1); ++i)
        {
            const ffe_t sum = AddModQ(GFLog[temp[i] ^ 1], temp[m]);
            temp[i] = mulE(temp[i], sum);
        }
    }

    for (unsigned i = 0; i < kFieldSize; ++i)
        skewVec[i] = GFLog[skewVec[i]];

#ifdef LEO_EXPERIMENT_EXTRA_MULS
    temp[0] = kModulus - temp[0];

    for (unsigned i = 1; i < (kGFBits - 1); ++i)
        temp[i] = (kModulus - temp[i] + temp[i - 1]) % kModulus;

    B[0] = 0;
    for (unsigned i = 0; i < (kGFBits - 1); ++i)
    {
        const unsigned depart = ((unsigned)1 << i);

        for (unsigned j = 0; j < depart; ++j)
            B[j + depart] = (B[j] + temp[i]) % kModulus;
    }
#endif

    for (unsigned i = 0; i < kFieldSize; ++i)
        log_walsh[i] = GFLog[i];

    log_walsh[0] = 0;

    FWHT(log_walsh, kGFBits);
}


//------------------------------------------------------------------------------
// Encoder

// Encoding alg for k/n<0.5: message is a power of two
static void encodeL(ffe_t* data, const unsigned k, ffe_t* codeword)
{
    memcpy(codeword, data, sizeof(ffe_t) * k);

    IFLT(codeword, k, 0);

    for (unsigned i = k; i < kFieldSize; i += k)
    {
        memcpy(&codeword[i], codeword, sizeof(ffe_t) * k);

        FLT(&codeword[i], k, i, k);
    }

    memcpy(codeword, data, sizeof(ffe_t) * k);
}

// Encoding alg for k/n>0.5: parity is a power of two.
// data: message array. parity: parity array. mem: buffer(size>= n-k)
static void encodeH(const ffe_t* data, const unsigned m, const unsigned original_count, ffe_t* parity, ffe_t* mem)
{
    // Note: Assumes data is padded with zeroes out to the next multiple of m

    memcpy(parity, data, m * sizeof(ffe_t));
    IFLT(parity, m, m);

    for (unsigned i = m; i < original_count; i += m)
    {
        memcpy(mem, data + i, m * sizeof(ffe_t));
        IFLT(mem, m, m + i);
        for (unsigned j = 0; j < m; ++j)
            parity[j] ^= mem[j];
    }

    FLT(parity, m, 0, m);
}


//------------------------------------------------------------------------------
// Decoder

static void decode(ffe_t* codeword, const unsigned m, const unsigned original_count, const unsigned n, const bool* erasure)
{
    fwht_t log_walsh2[kFieldSize];

    // Compute the evaluations of the error locator polynomial
    for (unsigned i = 0; i < kFieldSize; ++i)
        log_walsh2[i] = erasure[i] ? 1 : 0;

    FWHT(log_walsh2, kGFBits);

    for (unsigned i = 0; i < kFieldSize; ++i)
        log_walsh2[i] = ((unsigned)log_walsh2[i] * (unsigned)log_walsh[i]) % kModulus;

    FWHT(log_walsh2, kGFBits);

    // k2 can be replaced with k
    //const unsigned k2 = kFieldSize;
    //const unsigned k2 = k; // cannot actually be replaced with k.  maybe for encodeL() only?

    for (unsigned i = 0; i < m + original_count; ++i)
    {
        if (erasure[i])
        {
            codeword[i] = 0;
        }
        else
        {
            codeword[i] = mulE(codeword[i], log_walsh2[i]);
        }
    }
    for (unsigned i = m + original_count; i < n; ++i)
        codeword[i] = 0;

    IFLT(codeword, n, 0);

    // Note: This is not needed to recover successfully...
#ifdef LEO_EXPERIMENT_EXTRA_MULS
    // formal derivative
    // Note: Preserves zeroes on the right
    for (unsigned i = 0; i < m + original_count; i += 2)
    {
        codeword[i] = mulE(codeword[i], kModulus - B[i >> 1]);
        codeword[i + 1] = mulE(codeword[i + 1], kModulus - B[i >> 1]);
    }
#endif

    formal_derivative(codeword, n);

#ifdef LEO_EXPERIMENT_EXTRA_MULS
    // Note: Preserves zeroes on the right
    for (unsigned i = 0; i < m + original_count; i += 2)
    {
        codeword[i] = mulE(codeword[i], B[i >> 1]);
        codeword[i + 1] = mulE(codeword[i + 1], B[i >> 1]);
    }
#endif

    FLT(codeword, n, 0, m + original_count);

    for (unsigned i = 0; i < kFieldSize; ++i)
    {
        if (erasure[i])
        {
            codeword[i] = mulE(codeword[i], kModulus - log_walsh2[i]);
        }
    }
}


#ifdef _MSC_VER
#include <intrin.h>
#endif

// Returns highest bit index 0..63 where the first non-zero bit is found
// Precondition: x != 0
LEO_FORCE_INLINE unsigned LastNonzeroBit64(uint64_t x)
{
#ifdef _MSC_VER
#ifdef _WIN64
    unsigned long index;
    // Note: Ignoring result because x != 0
    _BitScanReverse64(&index, x);
    return (unsigned)index;
#else
    unsigned long index;
    if (0 != _BitScanReverse(&index, (uint32_t)x))
        return (unsigned)index;
    // Note: Ignoring result because x != 0
    _BitScanReverse(&index, (uint32_t)(x >> 32));
    return (unsigned)index + 32;
#endif
#else
    // Note: Ignoring return value of 0 because x != 0
    return 63 - (unsigned)__builtin_clzll(x);
#endif
}


//------------------------------------------------------------------------------
// Test Application

void test(unsigned original_count, unsigned recovery_count, unsigned seed)
{
    unsigned m = 2UL << LastNonzeroBit64(recovery_count - 1);
    unsigned n = 2UL << LastNonzeroBit64(m + original_count - 1);

    srand(seed);

    //-----------Generating message----------

    // Message array
    ffe_t data[kFieldSize] = {0};

    // Filled with random numbers
    for (unsigned i = m; i < m + original_count; ++i)
        data[i] = (ffe_t)rand();


    //---------encoding----------

    ffe_t codeword[kFieldSize] = {};
    // First m codewords are for the parity data
    encodeH(data + m, m, original_count, data, codeword);
    //encodeL(data, k, codeword); // does not seem to work with any input?  what else needs to change?

    memcpy(codeword, data, sizeof(ffe_t) * kFieldSize);


    //--------erasure simulation---------

    // Array indicating erasures
    bool erasure[kFieldSize] = {
        false
    };

    // Tag the first "recovery_count" elements as erasures
    for (unsigned i = m; i < m + recovery_count; ++i)
        erasure[i] = true;

#if 0
    // permuting the erasure array
    for (unsigned i = m + original_count - 1; i > 0; --i)
    {
        unsigned pos = rand() % (i + 1);

        if (i != pos)
        {
            bool tmp = erasure[i];
            erasure[i] = erasure[pos];
            erasure[pos] = tmp;
        }
    }
#endif


    //---------main processing----------
    decode(codeword, m, original_count, n, erasure);

    // Check the correctness of the result
    for (unsigned i = 0; i < kFieldSize; ++i)
    {
        if (erasure[i])
        {
            if (data[i] != codeword[i])
            {
                printf("Decoding Error with seed = %d!\n", seed);
                LEO_DEBUG_BREAK;
                return;
            }
        }
    }

    printf(":D ");
}


//------------------------------------------------------------------------------
// Entrypoint

int main(int argc, char **argv)
{
    // Fill GFLog table and GFExp table
    InitField();

    // Compute factors used in erasure decoder
    InitFieldOperations();

    unsigned seed = (unsigned)time(NULL);
    for (;;)
    {
#ifdef LEO_SHORT_FIELD
        const unsigned input_count = 32;
        const unsigned recovery_count = 16;
#else // LEO_SHORT_FIELD
        const unsigned input_count = 32768;
        const unsigned recovery_count = 32768;
#endif // LEO_SHORT_FIELD

        test(input_count, recovery_count, seed);

        ++seed;
    }

    return 0;
}
