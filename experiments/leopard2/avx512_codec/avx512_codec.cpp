/*
 * Standalone AVX-512BW/VBMI experiment for Leopard's legacy GF8 wire field.
 * This intentionally does not call production Leopard arithmetic.
 */
#include <immintrin.h>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

constexpr unsigned kPolynomial = 0x11d;
constexpr std::array<uint8_t, 8> kCantorBasis =
    {{1, 214, 152, 146, 86, 200, 88, 230}};
constexpr unsigned kSamples = 9;
constexpr uint64_t kTargetTraffic = 32ULL * 1024ULL * 1024ULL;

#if defined(__GNUC__) || defined(__clang__)
#define NOINLINE __attribute__((noinline))
#define AVX2_FN __attribute__((target("avx2"), noinline))
#define AVX512BW_FN __attribute__((target("avx512f,avx512bw"), noinline))
#define AVX512VBMI_FN \
    __attribute__((target("avx512f,avx512bw,avx512vbmi"), noinline))
#define AVX2_INLINE __attribute__((target("avx2"), always_inline)) inline
#define AVX512BW_INLINE \
    __attribute__((target("avx512f,avx512bw"), always_inline)) inline
#define AVX512VBMI_INLINE \
    __attribute__((target("avx512f,avx512bw,avx512vbmi"), always_inline)) inline
#else
#error This standalone experiment requires GCC or Clang target attributes
#endif

void Require(bool condition, const std::string& message)
{
    if (!condition)
        throw std::runtime_error(message);
}

uint8_t PolynomialMultiply(uint8_t a, uint8_t b)
{
    unsigned x = a;
    unsigned y = b;
    unsigned result = 0;
    while (y != 0)
    {
        if ((y & 1U) != 0)
            result ^= x;
        y >>= 1U;
        x <<= 1U;
        if ((x & 0x100U) != 0)
            x ^= kPolynomial;
    }
    return static_cast<uint8_t>(result);
}

uint8_t CantorToPolynomial(uint8_t value)
{
    uint8_t result = 0;
    for (unsigned bit = 0; bit < 8; ++bit)
        if ((value & static_cast<uint8_t>(1U << bit)) != 0)
            result ^= kCantorBasis[bit];
    return result;
}

struct Tables
{
    alignas(64) std::array<std::array<uint8_t, 256>, 256> multiply{};
    alignas(64) std::array<std::array<uint8_t, 64>, 256> nibble_low{};
    alignas(64) std::array<std::array<uint8_t, 64>, 256> nibble_high{};
};

Tables g_tables;

void BuildTables()
{
    std::array<uint8_t, 256> polynomial_to_cantor{};
    for (unsigned i = 0; i < 256; ++i)
        polynomial_to_cantor[CantorToPolynomial(static_cast<uint8_t>(i))] =
            static_cast<uint8_t>(i);

    for (unsigned multiplier = 0; multiplier < 256; ++multiplier)
    {
        for (unsigned input = 0; input < 256; ++input)
        {
            const uint8_t product = PolynomialMultiply(
                CantorToPolynomial(static_cast<uint8_t>(multiplier)),
                CantorToPolynomial(static_cast<uint8_t>(input)));
            g_tables.multiply[multiplier][input] =
                polynomial_to_cantor[product];
        }
        for (unsigned i = 0; i < 64; ++i)
        {
            const unsigned nibble = i & 15U;
            g_tables.nibble_low[multiplier][i] =
                g_tables.multiply[multiplier][nibble];
            g_tables.nibble_high[multiplier][i] =
                g_tables.multiply[multiplier][nibble << 4U];
        }
    }
}

uint64_t Fnv1a(const uint8_t* data, size_t bytes)
{
    uint64_t hash = 1469598103934665603ULL;
    for (size_t i = 0; i < bytes; ++i)
    {
        hash ^= data[i];
        hash *= 1099511628211ULL;
    }
    return hash;
}

NOINLINE void ScalarCopy(uint8_t* out, const uint8_t* in, size_t bytes)
{
    for (size_t i = 0; i < bytes; ++i)
        out[i] = in[i];
}

NOINLINE void ScalarXor(uint8_t* out, const uint8_t* a, const uint8_t* b,
                        size_t bytes)
{
    for (size_t i = 0; i < bytes; ++i)
        out[i] = static_cast<uint8_t>(a[i] ^ b[i]);
}

NOINLINE void ScalarMul(uint8_t* out, const uint8_t* in, size_t bytes,
                        uint8_t coefficient)
{
    const uint8_t* row = g_tables.multiply[coefficient].data();
    for (size_t i = 0; i < bytes; ++i)
        out[i] = row[in[i]];
}

inline void ScalarButterfly(uint8_t& x, uint8_t& y, uint8_t coefficient)
{
    y ^= x;
    x ^= g_tables.multiply[coefficient][y];
}

constexpr std::array<uint8_t, 7> kFactors = {{2, 53, 91, 197, 7, 144, 233}};

NOINLINE void ScalarRadix4(uint8_t* const* work, size_t bytes)
{
    for (size_t i = 0; i < bytes; ++i)
    {
        uint8_t x0 = work[0][i], x1 = work[1][i];
        uint8_t x2 = work[2][i], x3 = work[3][i];
        ScalarButterfly(x0, x1, kFactors[0]);
        ScalarButterfly(x2, x3, kFactors[1]);
        ScalarButterfly(x0, x2, kFactors[4]);
        ScalarButterfly(x1, x3, kFactors[4]);
        work[0][i] = x0; work[1][i] = x1;
        work[2][i] = x2; work[3][i] = x3;
    }
}

NOINLINE void ScalarRadix8(uint8_t* const* work, size_t bytes)
{
    for (size_t i = 0; i < bytes; ++i)
    {
        std::array<uint8_t, 8> x{};
        for (unsigned j = 0; j < 8; ++j) x[j] = work[j][i];
        ScalarButterfly(x[0], x[1], kFactors[0]);
        ScalarButterfly(x[2], x[3], kFactors[1]);
        ScalarButterfly(x[4], x[5], kFactors[2]);
        ScalarButterfly(x[6], x[7], kFactors[3]);
        ScalarButterfly(x[0], x[2], kFactors[4]);
        ScalarButterfly(x[1], x[3], kFactors[4]);
        ScalarButterfly(x[4], x[6], kFactors[5]);
        ScalarButterfly(x[5], x[7], kFactors[5]);
        for (unsigned j = 0; j < 4; ++j)
            ScalarButterfly(x[j], x[j + 4], kFactors[6]);
        for (unsigned j = 0; j < 8; ++j) work[j][i] = x[j];
    }
}

constexpr std::array<std::array<uint8_t, 8>, 4> kCodecCoefficients = {{
    {{1, 2, 3, 5, 8, 13, 21, 34}},
    {{55, 89, 144, 233, 1, 7, 0, 197}},
    {{11, 29, 47, 83, 131, 173, 211, 251}},
    {{0, 1, 254, 127, 63, 31, 15, 9}}
}};

NOINLINE void ScalarCodec(const uint8_t* const* input, uint8_t* const* parity,
                          size_t bytes)
{
    for (size_t offset = 0; offset < bytes; ++offset)
    {
        std::array<uint8_t, 4> p{};
        for (unsigned row = 0; row < 4; ++row)
            for (unsigned column = 0; column < 8; ++column)
                p[row] ^= g_tables.multiply[kCodecCoefficients[row][column]]
                                             [input[column][offset]];
        ScalarButterfly(p[0], p[1], kFactors[0]);
        ScalarButterfly(p[2], p[3], kFactors[1]);
        ScalarButterfly(p[0], p[2], kFactors[4]);
        ScalarButterfly(p[1], p[3], kFactors[4]);
        for (unsigned row = 0; row < 4; ++row)
            parity[row][offset] = p[row];
    }
}

AVX2_INLINE __m256i Mul256(__m256i x, uint8_t coefficient)
{
    if (coefficient == 0) return _mm256_setzero_si256();
    if (coefficient == 1) return x;
    const __m256i lo_table = _mm256_load_si256(
        reinterpret_cast<const __m256i*>(g_tables.nibble_low[coefficient].data()));
    const __m256i hi_table = _mm256_load_si256(
        reinterpret_cast<const __m256i*>(g_tables.nibble_high[coefficient].data()));
    const __m256i mask = _mm256_set1_epi8(0x0f);
    const __m256i lo = _mm256_shuffle_epi8(lo_table, _mm256_and_si256(x, mask));
    const __m256i hi_index = _mm256_and_si256(_mm256_srli_epi64(x, 4), mask);
    return _mm256_xor_si256(lo, _mm256_shuffle_epi8(hi_table, hi_index));
}

AVX512BW_INLINE __m512i Mul512BW(__m512i x, uint8_t coefficient)
{
    if (coefficient == 0) return _mm512_setzero_si512();
    if (coefficient == 1) return x;
    const __m512i lo_table = _mm512_load_si512(g_tables.nibble_low[coefficient].data());
    const __m512i hi_table = _mm512_load_si512(g_tables.nibble_high[coefficient].data());
    const __m512i mask = _mm512_set1_epi8(0x0f);
    const __m512i lo = _mm512_shuffle_epi8(lo_table, _mm512_and_si512(x, mask));
    const __m512i hi_index = _mm512_and_si512(_mm512_srli_epi64(x, 4), mask);
    return _mm512_xor_si512(lo, _mm512_shuffle_epi8(hi_table, hi_index));
}

AVX512VBMI_INLINE __m512i Mul512VBMI(__m512i x, uint8_t coefficient)
{
    if (coefficient == 0) return _mm512_setzero_si512();
    if (coefficient == 1) return x;
    const uint8_t* row = g_tables.multiply[coefficient].data();
    const __m512i t0 = _mm512_load_si512(row);
    const __m512i t1 = _mm512_load_si512(row + 64);
    const __m512i t2 = _mm512_load_si512(row + 128);
    const __m512i t3 = _mm512_load_si512(row + 192);
    const __m512i high = _mm512_and_si512(x, _mm512_set1_epi8(static_cast<char>(0xc0)));
    __m512i out = _mm512_permutexvar_epi8(x, t0);
    out = _mm512_mask_permutexvar_epi8(out,
        _mm512_cmpeq_epi8_mask(high, _mm512_set1_epi8(0x40)), x, t1);
    out = _mm512_mask_permutexvar_epi8(out,
        _mm512_cmpeq_epi8_mask(high, _mm512_set1_epi8(static_cast<char>(0x80))), x, t2);
    out = _mm512_mask_permutexvar_epi8(out,
        _mm512_cmpeq_epi8_mask(high, _mm512_set1_epi8(static_cast<char>(0xc0))), x, t3);
    return out;
}

AVX2_FN void Avx2Copy(uint8_t* out, const uint8_t* in, size_t bytes)
{
    size_t i = 0;
    for (; i + 32 <= bytes; i += 32)
        _mm256_storeu_si256(reinterpret_cast<__m256i*>(out + i),
            _mm256_loadu_si256(reinterpret_cast<const __m256i*>(in + i)));
    ScalarCopy(out + i, in + i, bytes - i);
}

AVX512BW_FN void Avx512Copy(uint8_t* out, const uint8_t* in, size_t bytes)
{
    size_t i = 0;
    for (; i + 64 <= bytes; i += 64)
        _mm512_storeu_si512(out + i, _mm512_loadu_si512(in + i));
    ScalarCopy(out + i, in + i, bytes - i);
}

AVX2_FN void Avx2Xor(uint8_t* out, const uint8_t* a, const uint8_t* b, size_t bytes)
{
    size_t i = 0;
    for (; i + 32 <= bytes; i += 32)
        _mm256_storeu_si256(reinterpret_cast<__m256i*>(out + i),
            _mm256_xor_si256(
                _mm256_loadu_si256(reinterpret_cast<const __m256i*>(a + i)),
                _mm256_loadu_si256(reinterpret_cast<const __m256i*>(b + i))));
    ScalarXor(out + i, a + i, b + i, bytes - i);
}

AVX512BW_FN void Avx512Xor(uint8_t* out, const uint8_t* a, const uint8_t* b,
                           size_t bytes)
{
    size_t i = 0;
    for (; i + 64 <= bytes; i += 64)
        _mm512_storeu_si512(out + i,
            _mm512_xor_si512(_mm512_loadu_si512(a + i), _mm512_loadu_si512(b + i)));
    ScalarXor(out + i, a + i, b + i, bytes - i);
}

AVX2_FN void Avx2Mul(uint8_t* out, const uint8_t* in, size_t bytes, uint8_t c)
{
    size_t i = 0;
    for (; i + 32 <= bytes; i += 32)
        _mm256_storeu_si256(reinterpret_cast<__m256i*>(out + i),
            Mul256(_mm256_loadu_si256(reinterpret_cast<const __m256i*>(in + i)), c));
    ScalarMul(out + i, in + i, bytes - i, c);
}

AVX512BW_FN void Avx512BWMul(uint8_t* out, const uint8_t* in, size_t bytes, uint8_t c)
{
    size_t i = 0;
    for (; i + 64 <= bytes; i += 64)
        _mm512_storeu_si512(out + i, Mul512BW(_mm512_loadu_si512(in + i), c));
    ScalarMul(out + i, in + i, bytes - i, c);
}

AVX512VBMI_FN void Avx512VBMIMul(uint8_t* out, const uint8_t* in, size_t bytes,
                                 uint8_t c)
{
    size_t i = 0;
    for (; i + 64 <= bytes; i += 64)
        _mm512_storeu_si512(out + i, Mul512VBMI(_mm512_loadu_si512(in + i), c));
    ScalarMul(out + i, in + i, bytes - i, c);
}

AVX2_FN void Avx2Radix4(uint8_t* const* w, size_t bytes)
{
    size_t i = 0;
    for (; i + 32 <= bytes; i += 32)
    {
        __m256i x0 = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(w[0] + i));
        __m256i x1 = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(w[1] + i));
        __m256i x2 = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(w[2] + i));
        __m256i x3 = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(w[3] + i));
        x1 ^= x0; x0 ^= Mul256(x1, kFactors[0]);
        x3 ^= x2; x2 ^= Mul256(x3, kFactors[1]);
        x2 ^= x0; x0 ^= Mul256(x2, kFactors[4]);
        x3 ^= x1; x1 ^= Mul256(x3, kFactors[4]);
        _mm256_storeu_si256(reinterpret_cast<__m256i*>(w[0] + i), x0);
        _mm256_storeu_si256(reinterpret_cast<__m256i*>(w[1] + i), x1);
        _mm256_storeu_si256(reinterpret_cast<__m256i*>(w[2] + i), x2);
        _mm256_storeu_si256(reinterpret_cast<__m256i*>(w[3] + i), x3);
    }
    std::array<uint8_t*, 4> tail = {{w[0]+i,w[1]+i,w[2]+i,w[3]+i}};
    ScalarRadix4(tail.data(), bytes - i);
}

AVX512BW_FN void Avx512BWRadix4(uint8_t* const* w, size_t bytes)
{
    size_t i=0;
    for (;i+64<=bytes;i+=64) {
        __m512i x0=_mm512_loadu_si512(w[0]+i),x1=_mm512_loadu_si512(w[1]+i);
        __m512i x2=_mm512_loadu_si512(w[2]+i),x3=_mm512_loadu_si512(w[3]+i);
        x1^=x0;x0^=Mul512BW(x1,kFactors[0]);x3^=x2;x2^=Mul512BW(x3,kFactors[1]);
        x2^=x0;x0^=Mul512BW(x2,kFactors[4]);x3^=x1;x1^=Mul512BW(x3,kFactors[4]);
        _mm512_storeu_si512(w[0]+i,x0);_mm512_storeu_si512(w[1]+i,x1);
        _mm512_storeu_si512(w[2]+i,x2);_mm512_storeu_si512(w[3]+i,x3);
    }
    std::array<uint8_t*,4> t={{w[0]+i,w[1]+i,w[2]+i,w[3]+i}};ScalarRadix4(t.data(),bytes-i);
}

AVX512VBMI_FN void Avx512VBMIRadix4(uint8_t* const* w, size_t bytes)
{
    size_t i=0;
    for (;i+64<=bytes;i+=64) {
        __m512i x0=_mm512_loadu_si512(w[0]+i),x1=_mm512_loadu_si512(w[1]+i);
        __m512i x2=_mm512_loadu_si512(w[2]+i),x3=_mm512_loadu_si512(w[3]+i);
        x1^=x0;x0^=Mul512VBMI(x1,kFactors[0]);x3^=x2;x2^=Mul512VBMI(x3,kFactors[1]);
        x2^=x0;x0^=Mul512VBMI(x2,kFactors[4]);x3^=x1;x1^=Mul512VBMI(x3,kFactors[4]);
        _mm512_storeu_si512(w[0]+i,x0);_mm512_storeu_si512(w[1]+i,x1);
        _mm512_storeu_si512(w[2]+i,x2);_mm512_storeu_si512(w[3]+i,x3);
    }
    std::array<uint8_t*,4> t={{w[0]+i,w[1]+i,w[2]+i,w[3]+i}};ScalarRadix4(t.data(),bytes-i);
}

#define RADIX8_BODY(VEC, WIDTH, LOAD, STORE, MUL) \
    size_t i=0; for(;i+WIDTH<=bytes;i+=WIDTH){ \
        VEC x[8]; for(unsigned j=0;j<8;++j)x[j]=LOAD(w[j]+i); \
        x[1]^=x[0];x[0]^=MUL(x[1],kFactors[0]);x[3]^=x[2];x[2]^=MUL(x[3],kFactors[1]); \
        x[5]^=x[4];x[4]^=MUL(x[5],kFactors[2]);x[7]^=x[6];x[6]^=MUL(x[7],kFactors[3]); \
        x[2]^=x[0];x[0]^=MUL(x[2],kFactors[4]);x[3]^=x[1];x[1]^=MUL(x[3],kFactors[4]); \
        x[6]^=x[4];x[4]^=MUL(x[6],kFactors[5]);x[7]^=x[5];x[5]^=MUL(x[7],kFactors[5]); \
        for(unsigned j=0;j<4;++j){x[j+4]^=x[j];x[j]^=MUL(x[j+4],kFactors[6]);} \
        for(unsigned j=0;j<8;++j)STORE(w[j]+i,x[j]);} \
    std::array<uint8_t*,8> t;for(unsigned j=0;j<8;++j)t[j]=w[j]+i;ScalarRadix8(t.data(),bytes-i)

AVX2_INLINE __m256i Load256(const uint8_t* p){return _mm256_loadu_si256(reinterpret_cast<const __m256i*>(p));}
AVX2_INLINE void Store256(uint8_t* p,__m256i x){_mm256_storeu_si256(reinterpret_cast<__m256i*>(p),x);}
AVX512BW_INLINE __m512i Load512BW(const uint8_t* p){return _mm512_loadu_si512(p);}
AVX512BW_INLINE void Store512BW(uint8_t* p,__m512i x){_mm512_storeu_si512(p,x);}
AVX512VBMI_INLINE __m512i Load512VBMI(const uint8_t* p){return _mm512_loadu_si512(p);}
AVX512VBMI_INLINE void Store512VBMI(uint8_t* p,__m512i x){_mm512_storeu_si512(p,x);}

AVX2_FN void Avx2Radix8(uint8_t* const* w,size_t bytes){RADIX8_BODY(__m256i,32,Load256,Store256,Mul256);}
AVX512BW_FN void Avx512BWRadix8(uint8_t* const* w,size_t bytes){RADIX8_BODY(__m512i,64,Load512BW,Store512BW,Mul512BW);}
AVX512VBMI_FN void Avx512VBMIRadix8(uint8_t* const* w,size_t bytes){RADIX8_BODY(__m512i,64,Load512VBMI,Store512VBMI,Mul512VBMI);}

#define CODEC_BODY(VEC, WIDTH, LOAD, STORE, MUL) \
    size_t i=0;for(;i+WIDTH<=bytes;i+=WIDTH){VEC p[4]; \
        for(unsigned r=0;r<4;++r){p[r]=MUL(LOAD(input[0]+i),kCodecCoefficients[r][0]); \
            for(unsigned c=1;c<8;++c)p[r]^=MUL(LOAD(input[c]+i),kCodecCoefficients[r][c]);} \
        p[1]^=p[0];p[0]^=MUL(p[1],kFactors[0]);p[3]^=p[2];p[2]^=MUL(p[3],kFactors[1]); \
        p[2]^=p[0];p[0]^=MUL(p[2],kFactors[4]);p[3]^=p[1];p[1]^=MUL(p[3],kFactors[4]); \
        for(unsigned r=0;r<4;++r)STORE(parity[r]+i,p[r]);} \
    std::array<const uint8_t*,8> ti;std::array<uint8_t*,4> tp; \
    for(unsigned j=0;j<8;++j){ti[j]=input[j]+i;} \
    for(unsigned j=0;j<4;++j){tp[j]=parity[j]+i;} \
    ScalarCodec(ti.data(),tp.data(),bytes-i)

AVX2_FN void Avx2Codec(const uint8_t* const* input,uint8_t* const* parity,size_t bytes){CODEC_BODY(__m256i,32,Load256,Store256,Mul256);}
AVX512BW_FN void Avx512BWCodec(const uint8_t* const* input,uint8_t* const* parity,size_t bytes){CODEC_BODY(__m512i,64,Load512BW,Store512BW,Mul512BW);}
AVX512VBMI_FN void Avx512VBMICodec(const uint8_t* const* input,uint8_t* const* parity,size_t bytes){CODEC_BODY(__m512i,64,Load512VBMI,Store512VBMI,Mul512VBMI);}

using CopyFn=void(*)(uint8_t*,const uint8_t*,size_t);
using XorFn=void(*)(uint8_t*,const uint8_t*,const uint8_t*,size_t);
using MulFn=void(*)(uint8_t*,const uint8_t*,size_t,uint8_t);
using RadixFn=void(*)(uint8_t* const*,size_t);
using CodecFn=void(*)(const uint8_t* const*,uint8_t* const*,size_t);

struct Backend { const char* name; CopyFn copy; XorFn xorf; MulFn mul; RadixFn r4; RadixFn r8; CodecFn codec; };

const std::array<Backend,4> kBackends = {{
    {"scalar",ScalarCopy,ScalarXor,ScalarMul,ScalarRadix4,ScalarRadix8,ScalarCodec},
    {"avx2",Avx2Copy,Avx2Xor,Avx2Mul,Avx2Radix4,Avx2Radix8,Avx2Codec},
    {"avx512bw",Avx512Copy,Avx512Xor,Avx512BWMul,Avx512BWRadix4,Avx512BWRadix8,Avx512BWCodec},
    {"avx512vbmi",Avx512Copy,Avx512Xor,Avx512VBMIMul,Avx512VBMIRadix4,Avx512VBMIRadix8,Avx512VBMICodec}
}};

uint64_t NextRandom(uint64_t& state){state^=state<<13;state^=state>>7;state^=state<<17;return state;}

void Fill(std::vector<uint8_t>& v,uint64_t& state){for(uint8_t& x:v)x=static_cast<uint8_t>(NextRandom(state));}

void Validate()
{
    Require(__builtin_cpu_supports("avx2"),"AVX2 unavailable");
    Require(__builtin_cpu_supports("avx512bw"),"AVX-512BW unavailable");
    Require(__builtin_cpu_supports("avx512vbmi"),"AVX-512VBMI unavailable");
    uint64_t products=0;
    std::array<uint8_t,256> inverse{};
    for(unsigned i=0;i<256;++i)inverse[CantorToPolynomial(static_cast<uint8_t>(i))]=static_cast<uint8_t>(i);
    for(unsigned a=0;a<256;++a)for(unsigned b=0;b<256;++b){
        const uint8_t expected=inverse[PolynomialMultiply(CantorToPolynomial(static_cast<uint8_t>(a)),CantorToPolynomial(static_cast<uint8_t>(b)))];
        Require(g_tables.multiply[a][b]==expected,"scalar table mismatch");++products;}

    std::vector<uint8_t> input(270),reference(270),candidate(270);
    for(unsigned i=0;i<270;++i)input[i]=static_cast<uint8_t>(i);
    uint64_t vector_products=0;
    for(unsigned c=0;c<256;++c){
        ScalarMul(reference.data()+3,input.data()+1,257,static_cast<uint8_t>(c));
        for(size_t b=1;b<kBackends.size();++b){
            std::fill(candidate.begin(),candidate.end(),0xa5);
            kBackends[b].mul(candidate.data()+5,input.data()+1,257,static_cast<uint8_t>(c));
            Require(std::memcmp(reference.data()+3,candidate.data()+5,257)==0,"SIMD multiplier mismatch");
            vector_products+=257;
        }}

    const std::array<size_t,22> sizes={{0,1,2,3,7,8,15,16,17,31,32,33,63,64,65,127,128,129,255,256,257,1023}};
    uint64_t tail_cases=0;uint64_t state=0x9182736455ULL;
    for(size_t bytes:sizes){
        std::array<std::vector<uint8_t>,8> base,ref,got;
        for(unsigned j=0;j<8;++j){base[j].resize(bytes+9);ref[j]=base[j];got[j]=base[j];Fill(base[j],state);ref[j]=base[j];got[j]=base[j];}
        for(size_t backend=1;backend<kBackends.size();++backend){
            std::vector<uint8_t> copy_ref(bytes+9),copy_got(bytes+9),xor_ref(bytes+9),xor_got(bytes+9);
            ScalarCopy(copy_ref.data()+2,base[0].data()+1,bytes);kBackends[backend].copy(copy_got.data()+3,base[0].data()+1,bytes);
            ScalarXor(xor_ref.data()+2,base[0].data()+1,base[1].data()+2,bytes);kBackends[backend].xorf(xor_got.data()+3,base[0].data()+1,base[1].data()+2,bytes);
            Require(std::memcmp(copy_ref.data()+2,copy_got.data()+3,bytes)==0,"copy tail mismatch");
            Require(std::memcmp(xor_ref.data()+2,xor_got.data()+3,bytes)==0,"xor tail mismatch");
            std::array<uint8_t*,8> rp{},gp{};for(unsigned j=0;j<8;++j){ref[j]=base[j];got[j]=base[j];rp[j]=ref[j].data()+1;gp[j]=got[j].data()+1;}
            ScalarRadix4(rp.data(),bytes);kBackends[backend].r4(gp.data(),bytes);
            for(unsigned j=0;j<4;++j)Require(std::memcmp(rp[j],gp[j],bytes)==0,"radix4 mismatch");
            for(unsigned j=0;j<8;++j){ref[j]=base[j];got[j]=base[j];rp[j]=ref[j].data()+1;gp[j]=got[j].data()+1;}
            ScalarRadix8(rp.data(),bytes);kBackends[backend].r8(gp.data(),bytes);
            for(unsigned j=0;j<8;++j)Require(std::memcmp(rp[j],gp[j],bytes)==0,"radix8 mismatch");
            std::array<const uint8_t*,8> in{};std::array<uint8_t*,4> ro{},go{};
            for(unsigned j=0;j<8;++j)in[j]=base[j].data()+1;
            for(unsigned j=0;j<4;++j){ref[j].assign(bytes+9,0);got[j].assign(bytes+9,0);ro[j]=ref[j].data()+2;go[j]=got[j].data()+3;}
            ScalarCodec(in.data(),ro.data(),bytes);kBackends[backend].codec(in.data(),go.data(),bytes);
            for(unsigned j=0;j<4;++j)Require(std::memcmp(ro[j],go[j],bytes)==0,"codec mismatch");
            ++tail_cases;
        }}
    const uint64_t hash=Fnv1a(g_tables.multiply[0].data(),256U*256U);
    std::cout<<"validation_passed,scalar_products="<<products<<",simd_products="<<vector_products
             <<",tail_backend_cases="<<tail_cases<<",table_fnv1a=0x"<<std::hex<<hash<<std::dec<<"\n";
}

struct Arena
{
    uint8_t* data=nullptr;size_t stride=0;unsigned count=0;
    Arena(unsigned n,size_t bytes):stride((bytes+63U)&~size_t(63U)),count(n){
        void* p=nullptr;Require(posix_memalign(&p,64,stride*count)==0,"posix_memalign failed");data=static_cast<uint8_t*>(p);}
    ~Arena(){std::free(data);}Arena(const Arena&)=delete;Arena& operator=(const Arena&)=delete;
    uint8_t* operator[](unsigned i){return data+stride*i;}
};

double ReadFrequencyMHz(int cpu)
{
    std::ifstream in("/sys/devices/system/cpu/cpu"+std::to_string(cpu)+"/cpufreq/scaling_cur_freq");
    double khz=0;in>>khz;return khz/1000.0;
}

double Median(std::vector<double> values){std::sort(values.begin(),values.end());return values[values.size()/2];}

template<class Function>
void EmitCell(const std::string& operation,const std::string& variant,size_t bytes,
              unsigned traffic_multiplier,int cpu,Function&& function,uint64_t& sink)
{
    const uint64_t logical=std::max<uint64_t>(1,static_cast<uint64_t>(bytes)*traffic_multiplier);
    uint64_t iterations=(kTargetTraffic+logical-1)/logical;iterations=std::max<uint64_t>(1,std::min<uint64_t>(iterations,1048576));
    function();function();
    std::vector<double> ns,cycles,before,after;ns.reserve(kSamples);cycles.reserve(kSamples);
    for(unsigned sample=0;sample<kSamples;++sample){
        before.push_back(ReadFrequencyMHz(cpu));unsigned aux=0;_mm_lfence();const uint64_t c0=__rdtsc();
        const auto t0=std::chrono::steady_clock::now();for(uint64_t j=0;j<iterations;++j)function();
        const auto t1=std::chrono::steady_clock::now();const uint64_t c1=__rdtscp(&aux);_mm_lfence();after.push_back(ReadFrequencyMHz(cpu));
        const double elapsed=std::chrono::duration<double,std::nano>(t1-t0).count();ns.push_back(elapsed/static_cast<double>(iterations));cycles.push_back(static_cast<double>(c1-c0)/static_cast<double>(iterations));}
    const double med=Median(ns);std::vector<double> deviations;for(double x:ns)deviations.push_back(std::abs(x-med));
    std::cout<<operation<<','<<variant<<','<<bytes<<','<<iterations<<','<<kSamples<<','<<std::fixed<<std::setprecision(3)
             <<med<<','<<Median(deviations)<<','<<*std::min_element(ns.begin(),ns.end())<<','<<*std::max_element(ns.begin(),ns.end())<<','
             <<Median(cycles)<<','<<(bytes==0?0.0:Median(cycles)/static_cast<double>(bytes))<<','
             <<(static_cast<double>(bytes)/med)<<','<<(static_cast<double>(logical)/med)<<','<<Median(before)<<','<<Median(after)<<','<<sink<<'\n';
}

void Benchmark(int cpu)
{
    const std::array<size_t,10> sizes={{64,256,1024,4096,16384,65536,262144,1048576,4194304,16777216}};
    const size_t maximum=sizes.back();Arena arena(16,maximum);uint64_t state=0x123456789abcdefULL;
    for(unsigned j=0;j<16;++j)for(size_t i=0;i<maximum;++i)arena[j][i]=static_cast<uint8_t>(NextRandom(state));
    std::cout<<"operation,variant,bytes,iterations,samples,median_ns,mad_ns,min_ns,max_ns,median_cycles,cycles_per_byte,shard_gb_s,logical_traffic_gb_s,freq_before_mhz,freq_after_mhz,checksum\n";
    uint64_t sink=0;
    for(size_t bytes:sizes)for(const Backend& backend:kBackends){
        EmitCell("copy",backend.name,bytes,2,cpu,[&]{backend.copy(arena[8],arena[0],bytes);},sink);sink^=arena[8][bytes-1];
        EmitCell("xor",backend.name,bytes,3,cpu,[&]{backend.xorf(arena[8],arena[0],arena[1],bytes);},sink);sink^=arena[8][bytes-1];
        EmitCell("mul",backend.name,bytes,2,cpu,[&]{backend.mul(arena[8],arena[0],bytes,173);},sink);sink^=arena[8][bytes-1];
        std::array<uint8_t*,8> work{};for(unsigned j=0;j<8;++j)work[j]=arena[8+j];
        EmitCell("radix4",backend.name,bytes,8,cpu,[&]{backend.r4(work.data(),bytes);},sink);sink^=arena[8][bytes-1];
        EmitCell("radix8",backend.name,bytes,16,cpu,[&]{backend.r8(work.data(),bytes);},sink);sink^=arena[15][bytes-1];
        std::array<const uint8_t*,8> input{};std::array<uint8_t*,4> parity{};for(unsigned j=0;j<8;++j)input[j]=arena[j];for(unsigned j=0;j<4;++j)parity[j]=arena[8+j];
        EmitCell("codec",backend.name,bytes,36,cpu,[&]{backend.codec(input.data(),parity.data(),bytes);},sink);sink^=arena[11][bytes-1];
    }
}

int ParseCpu(int argc,char** argv){for(int i=1;i+1<argc;++i)if(std::string(argv[i])=="--cpu")return std::stoi(argv[i+1]);return 0;}

} // namespace

int main(int argc,char** argv)
{
    try{
        BuildTables();const int cpu=ParseCpu(argc,argv);
        if(argc>=2&&std::string(argv[1])=="--validate"){Validate();return 0;}
        if(argc>=2&&std::string(argv[1])=="--benchmark"){Benchmark(cpu);return 0;}
        std::cerr<<"usage: "<<argv[0]<<" --validate | --benchmark --cpu N\n";return 2;
    }catch(const std::exception& e){std::cerr<<"error: "<<e.what()<<'\n';return 1;}
}
