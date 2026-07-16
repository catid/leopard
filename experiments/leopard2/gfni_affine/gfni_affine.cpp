/*
    Standalone Leopard2 GFNI-affine experiment.

    This file intentionally does not include or modify production Leopard
    internals.  It reconstructs the documented legacy GF8 representation and
    uses that independent scalar implementation as its correctness oracle.
*/

#include <immintrin.h>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#if defined(__linux__)
#include <sched.h>
#endif

namespace {

constexpr unsigned kOrder = 256;
constexpr unsigned kPolynomial = 0x11d;
constexpr uint8_t kCantorBasis[8] = {
    1, 214, 152, 146, 86, 200, 88, 230
};
constexpr unsigned kSamples = 9;

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
        if (y & 1u)
            result ^= x;
        y >>= 1;
        x <<= 1;
        if (x & 0x100u)
            x ^= kPolynomial;
    }
    return static_cast<uint8_t>(result);
}

uint8_t CantorToPolynomial(uint8_t value)
{
    uint8_t result = 0;
    for (unsigned bit = 0; bit < 8; ++bit)
        if (value & (1u << bit))
            result ^= kCantorBasis[bit];
    return result;
}

struct Tables
{
    std::array<uint8_t, 256> polynomial_to_cantor{};
    std::array<std::array<uint8_t, 256>, 256> multiply{};
    std::array<std::array<uint8_t, 16>, 256> nibble_low{};
    std::array<std::array<uint8_t, 16>, 256> nibble_high{};
    std::array<uint64_t, 256> affine_matrix{};
    std::array<std::array<unsigned, 8>, 8> operand_bit{};
    uint64_t to_aes_matrix = 0;
    uint64_t from_aes_matrix = 0;
    std::array<uint8_t, 256> aes_multiplier{};
    uint8_t aes_root = 0;
};

void BuildScalarTables(Tables& tables)
{
    for (unsigned cantor = 0; cantor < 256; ++cantor)
        tables.polynomial_to_cantor[CantorToPolynomial(
            static_cast<uint8_t>(cantor))] = static_cast<uint8_t>(cantor);

    for (unsigned multiplier = 0; multiplier < 256; ++multiplier)
        for (unsigned input = 0; input < 256; ++input)
        {
            const uint8_t product = PolynomialMultiply(
                CantorToPolynomial(static_cast<uint8_t>(multiplier)),
                CantorToPolynomial(static_cast<uint8_t>(input)));
            tables.multiply[multiplier][input] =
                tables.polynomial_to_cantor[product];
        }

    for (unsigned multiplier = 0; multiplier < 256; ++multiplier)
        for (unsigned nibble = 0; nibble < 16; ++nibble)
        {
            tables.nibble_low[multiplier][nibble] =
                tables.multiply[multiplier][nibble];
            tables.nibble_high[multiplier][nibble] =
                tables.multiply[multiplier][nibble << 4];
        }
}

__attribute__((target("gfni"), noinline))
uint8_t AffineProbe(uint8_t input, uint64_t matrix)
{
    const __m128i x = _mm_set1_epi8(static_cast<char>(input));
    const __m128i a = _mm_set1_epi64x(static_cast<long long>(matrix));
    const __m128i y = _mm_gf2p8affine_epi64_epi8(x, a, 0);
    return static_cast<uint8_t>(_mm_cvtsi128_si32(y));
}

void DiscoverMatrixBitOrder(Tables& tables)
{
    for (unsigned output_bit = 0; output_bit < 8; ++output_bit)
        for (unsigned input_bit = 0; input_bit < 8; ++input_bit)
            tables.operand_bit[output_bit][input_bit] = 64;

    for (unsigned operand_bit = 0; operand_bit < 64; ++operand_bit)
    {
        unsigned hits = 0;
        for (unsigned input_bit = 0; input_bit < 8; ++input_bit)
        {
            const uint8_t result = AffineProbe(
                static_cast<uint8_t>(1u << input_bit),
                UINT64_C(1) << operand_bit);
            if (result == 0)
                continue;
            Require((result & (result - 1)) == 0,
                "one matrix bit affected multiple output bits");
            const unsigned output_bit =
                static_cast<unsigned>(__builtin_ctz(result));
            Require(tables.operand_bit[output_bit][input_bit] == 64,
                "matrix bit-order probe produced a duplicate mapping");
            tables.operand_bit[output_bit][input_bit] = operand_bit;
            ++hits;
        }
        Require(hits == 1,
            "matrix bit-order probe did not produce one input/output pair");
    }

    for (unsigned output_bit = 0; output_bit < 8; ++output_bit)
        for (unsigned input_bit = 0; input_bit < 8; ++input_bit)
            Require(tables.operand_bit[output_bit][input_bit] < 64,
                "incomplete matrix bit-order mapping");
}

uint64_t EncodeLinearMap(
    const Tables& tables,
    const std::array<uint8_t, 8>& columns)
{
    uint64_t encoded = 0;
    for (unsigned input_bit = 0; input_bit < 8; ++input_bit)
        for (unsigned output_bit = 0; output_bit < 8; ++output_bit)
            if (columns[input_bit] & (1u << output_bit))
                encoded |= UINT64_C(1) <<
                    tables.operand_bit[output_bit][input_bit];
    return encoded;
}

void BuildAffineMatrices(Tables& tables)
{
    for (unsigned multiplier = 0; multiplier < 256; ++multiplier)
    {
        std::array<uint8_t, 8> columns{};
        for (unsigned input_bit = 0; input_bit < 8; ++input_bit)
            columns[input_bit] = tables.multiply[multiplier][1u << input_bit];
        tables.affine_matrix[multiplier] = EncodeLinearMap(tables, columns);
    }
}

uint8_t AESMultiply(uint8_t a, uint8_t b)
{
    unsigned x = a;
    unsigned y = b;
    unsigned result = 0;
    while (y != 0)
    {
        if (y & 1u)
            result ^= x;
        y >>= 1;
        x <<= 1;
        if (x & 0x100u)
            x ^= 0x11bu;
    }
    return static_cast<uint8_t>(result);
}

uint8_t AESPower(uint8_t value, unsigned exponent)
{
    uint8_t result = 1;
    while (exponent != 0)
    {
        if (exponent & 1u)
            result = AESMultiply(result, value);
        value = AESMultiply(value, value);
        exponent >>= 1;
    }
    return result;
}

void BuildAESIsomorphism(Tables& tables)
{
    // Find a root in the AES 0x11b field of Leopard's polynomial 0x11d:
    // x^8 + x^4 + x^3 + x^2 + 1.
    for (unsigned candidate = 1; candidate < 256; ++candidate)
    {
        const uint8_t x = static_cast<uint8_t>(candidate);
        const uint8_t evaluation = static_cast<uint8_t>(
            AESPower(x, 8) ^ AESPower(x, 4) ^ AESPower(x, 3) ^
            AESPower(x, 2) ^ 1);
        if (evaluation == 0)
        {
            tables.aes_root = x;
            break;
        }
    }
    Require(tables.aes_root != 0, "could not find AES-field root of 0x11d");

    std::array<uint8_t, 256> cantor_to_aes{};
    std::array<uint8_t, 256> aes_to_cantor{};
    std::array<uint8_t, 8> to_columns{};
    std::array<uint8_t, 8> from_columns{};
    for (unsigned cantor = 0; cantor < 256; ++cantor)
    {
        const uint8_t polynomial = CantorToPolynomial(static_cast<uint8_t>(cantor));
        uint8_t mapped = 0;
        for (unsigned bit = 0; bit < 8; ++bit)
            if (polynomial & (1u << bit))
                mapped ^= AESPower(tables.aes_root, bit);
        cantor_to_aes[cantor] = mapped;
        aes_to_cantor[mapped] = static_cast<uint8_t>(cantor);
    }
    for (unsigned bit = 0; bit < 8; ++bit)
    {
        to_columns[bit] = cantor_to_aes[1u << bit];
        from_columns[bit] = aes_to_cantor[1u << bit];
    }
    tables.to_aes_matrix = EncodeLinearMap(tables, to_columns);
    tables.from_aes_matrix = EncodeLinearMap(tables, from_columns);
    for (unsigned multiplier = 0; multiplier < 256; ++multiplier)
        tables.aes_multiplier[multiplier] = cantor_to_aes[multiplier];

    for (unsigned a = 0; a < 256; ++a)
        for (unsigned b = 0; b < 256; ++b)
            Require(aes_to_cantor[AESMultiply(
                cantor_to_aes[a], cantor_to_aes[b])] == tables.multiply[a][b],
                "AES-field basis isomorphism does not preserve multiplication");
}

uint64_t FNV1a64(const void* data, size_t bytes)
{
    const uint8_t* input = static_cast<const uint8_t*>(data);
    uint64_t hash = UINT64_C(1469598103934665603);
    for (size_t i = 0; i < bytes; ++i)
    {
        hash ^= input[i];
        hash *= UINT64_C(1099511628211);
    }
    return hash;
}

__attribute__((target("gfni"), noinline))
void RawGFNIXMM(uint8_t* output, const uint8_t* input, size_t bytes, uint64_t matrix)
{
    const __m128i a = _mm_set1_epi64x(static_cast<long long>(matrix));
    size_t offset = 0;
    for (; offset + 16 <= bytes; offset += 16)
    {
        const __m128i x = _mm_loadu_si128(
            reinterpret_cast<const __m128i*>(input + offset));
        const __m128i y = _mm_gf2p8affine_epi64_epi8(x, a, 0);
        _mm_storeu_si128(reinterpret_cast<__m128i*>(output + offset), y);
    }
    for (; offset < bytes; ++offset)
        output[offset] = AffineProbe(input[offset], matrix);
}

__attribute__((target("avx2,gfni"), noinline))
void RawGFNIYMM(uint8_t* output, const uint8_t* input, size_t bytes, uint64_t matrix)
{
    const __m256i a = _mm256_set1_epi64x(static_cast<long long>(matrix));
    size_t offset = 0;
    for (; offset + 32 <= bytes; offset += 32)
    {
        const __m256i x = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(input + offset));
        const __m256i y = _mm256_gf2p8affine_epi64_epi8(x, a, 0);
        _mm256_storeu_si256(reinterpret_cast<__m256i*>(output + offset), y);
    }
    for (; offset < bytes; ++offset)
        output[offset] = AffineProbe(input[offset], matrix);
}

__attribute__((target("avx512f,avx512bw,gfni"), noinline))
void RawGFNIZMM(uint8_t* output, const uint8_t* input, size_t bytes, uint64_t matrix)
{
    const __m512i a = _mm512_set1_epi64(static_cast<long long>(matrix));
    size_t offset = 0;
    for (; offset + 64 <= bytes; offset += 64)
    {
        const __m512i x = _mm512_loadu_si512(input + offset);
        const __m512i y = _mm512_gf2p8affine_epi64_epi8(x, a, 0);
        _mm512_storeu_si512(output + offset, y);
    }
    for (; offset < bytes; ++offset)
        output[offset] = AffineProbe(input[offset], matrix);
}

__attribute__((target("gfni"), noinline))
void RawGFNIMulXMM(uint8_t* output, const uint8_t* input, size_t bytes,
                   uint64_t to_aes, uint64_t from_aes, uint8_t multiplier)
{
    const __m128i to = _mm_set1_epi64x(static_cast<long long>(to_aes));
    const __m128i from = _mm_set1_epi64x(static_cast<long long>(from_aes));
    const __m128i m = _mm_set1_epi8(static_cast<char>(multiplier));
    size_t offset = 0;
    for (; offset + 16 <= bytes; offset += 16)
    {
        __m128i x = _mm_loadu_si128(reinterpret_cast<const __m128i*>(input + offset));
        x = _mm_gf2p8affine_epi64_epi8(x, to, 0);
        x = _mm_gf2p8mul_epi8(x, m);
        x = _mm_gf2p8affine_epi64_epi8(x, from, 0);
        _mm_storeu_si128(reinterpret_cast<__m128i*>(output + offset), x);
    }
}

__attribute__((target("avx2,gfni"), noinline))
void RawGFNIMulYMM(uint8_t* output, const uint8_t* input, size_t bytes,
                   uint64_t to_aes, uint64_t from_aes, uint8_t multiplier)
{
    const __m256i to = _mm256_set1_epi64x(static_cast<long long>(to_aes));
    const __m256i from = _mm256_set1_epi64x(static_cast<long long>(from_aes));
    const __m256i m = _mm256_set1_epi8(static_cast<char>(multiplier));
    size_t offset = 0;
    for (; offset + 32 <= bytes; offset += 32)
    {
        __m256i x = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(input + offset));
        x = _mm256_gf2p8affine_epi64_epi8(x, to, 0);
        x = _mm256_gf2p8mul_epi8(x, m);
        x = _mm256_gf2p8affine_epi64_epi8(x, from, 0);
        _mm256_storeu_si256(reinterpret_cast<__m256i*>(output + offset), x);
    }
}

__attribute__((target("avx512f,avx512bw,gfni"), noinline))
void RawGFNIMulZMM(uint8_t* output, const uint8_t* input, size_t bytes,
                   uint64_t to_aes, uint64_t from_aes, uint8_t multiplier)
{
    const __m512i to = _mm512_set1_epi64(static_cast<long long>(to_aes));
    const __m512i from = _mm512_set1_epi64(static_cast<long long>(from_aes));
    const __m512i m = _mm512_set1_epi8(static_cast<char>(multiplier));
    size_t offset = 0;
    for (; offset + 64 <= bytes; offset += 64)
    {
        __m512i x = _mm512_loadu_si512(input + offset);
        x = _mm512_gf2p8affine_epi64_epi8(x, to, 0);
        x = _mm512_gf2p8mul_epi8(x, m);
        x = _mm512_gf2p8affine_epi64_epi8(x, from, 0);
        _mm512_storeu_si512(output + offset, x);
    }
}

__attribute__((optimize("no-tree-vectorize"), noinline))
void ScalarTable(uint8_t* output, const uint8_t* input, size_t bytes,
                 const Tables& tables, uint8_t multiplier)
{
    if (multiplier == 0)
    {
        std::memset(output, 0, bytes);
        return;
    }
    if (multiplier == 1)
    {
        std::memcpy(output, input, bytes);
        return;
    }
    const uint8_t* table = tables.multiply[multiplier].data();
    for (size_t i = 0; i < bytes; ++i)
        output[i] = table[input[i]];
}

__attribute__((target("avx2"), noinline))
void AVX2Nibble(uint8_t* output, const uint8_t* input, size_t bytes,
                const Tables& tables, uint8_t multiplier)
{
    if (multiplier == 0)
    {
        std::memset(output, 0, bytes);
        return;
    }
    if (multiplier == 1)
    {
        std::memcpy(output, input, bytes);
        return;
    }
    const __m128i low128 = _mm_loadu_si128(reinterpret_cast<const __m128i*>(
        tables.nibble_low[multiplier].data()));
    const __m128i high128 = _mm_loadu_si128(reinterpret_cast<const __m128i*>(
        tables.nibble_high[multiplier].data()));
    const __m256i low_table = _mm256_broadcastsi128_si256(low128);
    const __m256i high_table = _mm256_broadcastsi128_si256(high128);
    const __m256i mask = _mm256_set1_epi8(15);
    size_t offset = 0;
    for (; offset + 32 <= bytes; offset += 32)
    {
        const __m256i x = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(input + offset));
        const __m256i lo = _mm256_and_si256(x, mask);
        const __m256i hi = _mm256_and_si256(_mm256_srli_epi64(x, 4), mask);
        const __m256i y = _mm256_xor_si256(
            _mm256_shuffle_epi8(low_table, lo),
            _mm256_shuffle_epi8(high_table, hi));
        _mm256_storeu_si256(reinterpret_cast<__m256i*>(output + offset), y);
    }
    for (; offset < bytes; ++offset)
        output[offset] = tables.multiply[multiplier][input[offset]];
}

__attribute__((target("avx512f,avx512bw"), noinline))
void AVX512Nibble(uint8_t* output, const uint8_t* input, size_t bytes,
                  const Tables& tables, uint8_t multiplier)
{
    if (multiplier == 0)
    {
        std::memset(output, 0, bytes);
        return;
    }
    if (multiplier == 1)
    {
        std::memcpy(output, input, bytes);
        return;
    }
    const __m128i low128 = _mm_loadu_si128(reinterpret_cast<const __m128i*>(
        tables.nibble_low[multiplier].data()));
    const __m128i high128 = _mm_loadu_si128(reinterpret_cast<const __m128i*>(
        tables.nibble_high[multiplier].data()));
    const __m512i low_table = _mm512_broadcast_i32x4(low128);
    const __m512i high_table = _mm512_broadcast_i32x4(high128);
    const __m512i mask = _mm512_set1_epi8(15);
    size_t offset = 0;
    for (; offset + 64 <= bytes; offset += 64)
    {
        const __m512i x = _mm512_loadu_si512(input + offset);
        const __m512i lo = _mm512_and_si512(x, mask);
        const __m512i hi = _mm512_and_si512(_mm512_srli_epi64(x, 4), mask);
        const __m512i y = _mm512_xor_si512(
            _mm512_shuffle_epi8(low_table, lo),
            _mm512_shuffle_epi8(high_table, hi));
        _mm512_storeu_si512(output + offset, y);
    }
    for (; offset < bytes; ++offset)
        output[offset] = tables.multiply[multiplier][input[offset]];
}

void GFNIXMM(uint8_t* output, const uint8_t* input, size_t bytes,
             const Tables& tables, uint8_t multiplier)
{
    if (multiplier == 0)
        std::memset(output, 0, bytes);
    else if (multiplier == 1)
        std::memcpy(output, input, bytes);
    else
        RawGFNIXMM(output, input, bytes, tables.affine_matrix[multiplier]);
}

void GFNIYMM(uint8_t* output, const uint8_t* input, size_t bytes,
             const Tables& tables, uint8_t multiplier)
{
    if (multiplier == 0)
        std::memset(output, 0, bytes);
    else if (multiplier == 1)
        std::memcpy(output, input, bytes);
    else
        RawGFNIYMM(output, input, bytes, tables.affine_matrix[multiplier]);
}

void GFNIZMM(uint8_t* output, const uint8_t* input, size_t bytes,
             const Tables& tables, uint8_t multiplier)
{
    if (multiplier == 0)
        std::memset(output, 0, bytes);
    else if (multiplier == 1)
        std::memcpy(output, input, bytes);
    else
        RawGFNIZMM(output, input, bytes, tables.affine_matrix[multiplier]);
}

void GFNIMulXMM(uint8_t* output, const uint8_t* input, size_t bytes,
                const Tables& tables, uint8_t multiplier)
{
    if (multiplier == 0)
        std::memset(output, 0, bytes);
    else if (multiplier == 1)
        std::memcpy(output, input, bytes);
    else
    {
        const size_t bulk = bytes & ~size_t(15);
        RawGFNIMulXMM(output, input, bulk, tables.to_aes_matrix,
            tables.from_aes_matrix, tables.aes_multiplier[multiplier]);
        for (size_t i = bulk; i < bytes; ++i)
            output[i] = tables.multiply[multiplier][input[i]];
    }
}

void GFNIMulYMM(uint8_t* output, const uint8_t* input, size_t bytes,
                const Tables& tables, uint8_t multiplier)
{
    if (multiplier == 0)
        std::memset(output, 0, bytes);
    else if (multiplier == 1)
        std::memcpy(output, input, bytes);
    else
    {
        const size_t bulk = bytes & ~size_t(31);
        RawGFNIMulYMM(output, input, bulk, tables.to_aes_matrix,
            tables.from_aes_matrix, tables.aes_multiplier[multiplier]);
        for (size_t i = bulk; i < bytes; ++i)
            output[i] = tables.multiply[multiplier][input[i]];
    }
}

void GFNIMulZMM(uint8_t* output, const uint8_t* input, size_t bytes,
                const Tables& tables, uint8_t multiplier)
{
    if (multiplier == 0)
        std::memset(output, 0, bytes);
    else if (multiplier == 1)
        std::memcpy(output, input, bytes);
    else
    {
        const size_t bulk = bytes & ~size_t(63);
        RawGFNIMulZMM(output, input, bulk, tables.to_aes_matrix,
            tables.from_aes_matrix, tables.aes_multiplier[multiplier]);
        for (size_t i = bulk; i < bytes; ++i)
            output[i] = tables.multiply[multiplier][input[i]];
    }
}

using Kernel = void (*)(uint8_t*, const uint8_t*, size_t,
                        const Tables&, uint8_t);

struct NamedKernel
{
    const char* name;
    Kernel function;
};

std::vector<NamedKernel> Kernels()
{
    return {
        { "scalar_table", ScalarTable },
        { "avx2_nibble", AVX2Nibble },
        { "avx512_nibble", AVX512Nibble },
        { "gfni_xmm", GFNIXMM },
        { "gfni_ymm", GFNIYMM },
        { "gfni_zmm", GFNIZMM },
        { "gfni_mul_xmm", GFNIMulXMM },
        { "gfni_mul_ymm", GFNIMulYMM },
        { "gfni_mul_zmm", GFNIMulZMM }
    };
}

void ValidateAll(const Tables& tables)
{
    std::array<uint8_t, 256> input{};
    std::array<uint8_t, 256> expected{};
    std::array<uint8_t, 256> output{};
    for (unsigned i = 0; i < 256; ++i)
        input[i] = static_cast<uint8_t>(i);

    for (unsigned multiplier = 0; multiplier < 256; ++multiplier)
    {
        for (unsigned input_value = 0; input_value < 256; ++input_value)
            expected[input_value] = tables.multiply[multiplier][input_value];

        RawGFNIXMM(output.data(), input.data(), input.size(),
            tables.affine_matrix[multiplier]);
        Require(output == expected, "XMM GFNI exhaustive mismatch");
        RawGFNIYMM(output.data(), input.data(), input.size(),
            tables.affine_matrix[multiplier]);
        Require(output == expected, "YMM GFNI exhaustive mismatch");
        RawGFNIZMM(output.data(), input.data(), input.size(),
            tables.affine_matrix[multiplier]);
        Require(output == expected, "ZMM GFNI exhaustive mismatch");
        RawGFNIMulXMM(output.data(), input.data(), input.size(),
            tables.to_aes_matrix, tables.from_aes_matrix,
            tables.aes_multiplier[multiplier]);
        Require(output == expected, "XMM GFNI multiply/isomorphism mismatch");
        RawGFNIMulYMM(output.data(), input.data(), input.size(),
            tables.to_aes_matrix, tables.from_aes_matrix,
            tables.aes_multiplier[multiplier]);
        Require(output == expected, "YMM GFNI multiply/isomorphism mismatch");
        RawGFNIMulZMM(output.data(), input.data(), input.size(),
            tables.to_aes_matrix, tables.from_aes_matrix,
            tables.aes_multiplier[multiplier]);
        Require(output == expected, "ZMM GFNI multiply/isomorphism mismatch");
        AVX2Nibble(output.data(), input.data(), input.size(), tables,
            static_cast<uint8_t>(multiplier));
        Require(output == expected, "AVX2 nibble exhaustive mismatch");
        AVX512Nibble(output.data(), input.data(), input.size(), tables,
            static_cast<uint8_t>(multiplier));
        Require(output == expected, "AVX512 nibble exhaustive mismatch");
    }

    std::array<uint8_t, 257> tail_input{};
    std::array<uint8_t, 257> tail_output{};
    for (unsigned i = 0; i < tail_input.size(); ++i)
        tail_input[i] = static_cast<uint8_t>(i * 73u + 19u);
    for (const NamedKernel& kernel : Kernels())
        for (uint8_t multiplier : { uint8_t(0), uint8_t(1), uint8_t(0x5d) })
        {
            kernel.function(tail_output.data(), tail_input.data(),
                tail_input.size(), tables, multiplier);
            for (unsigned i = 0; i < tail_input.size(); ++i)
                Require(tail_output[i] == tables.multiply[multiplier][tail_input[i]],
                    std::string(kernel.name) + " tail/specialization mismatch");
        }
}

void* AlignedAllocate(size_t bytes)
{
    void* pointer = nullptr;
    if (posix_memalign(&pointer, 64, bytes) != 0)
        throw std::bad_alloc();
    return pointer;
}

uint64_t ReadTSCStart()
{
    _mm_lfence();
    return __rdtsc();
}

uint64_t ReadTSCStop()
{
    unsigned auxiliary = 0;
    const uint64_t value = __rdtscp(&auxiliary);
    _mm_lfence();
    return value;
}

double Median(std::vector<double> values)
{
    std::sort(values.begin(), values.end());
    const size_t middle = values.size() / 2;
    return values.size() & 1u ? values[middle]
        : (values[middle - 1] + values[middle]) * 0.5;
}

struct Measurement
{
    double median_ns;
    double mad_ns;
    double median_cycles;
    size_t repetitions;
};

volatile uint64_t gSink = 0;

Measurement MeasureMultiply(
    const NamedKernel& kernel,
    const Tables& tables,
    uint8_t* output,
    const uint8_t* input,
    size_t bytes,
    uint8_t multiplier,
    unsigned passes)
{
    const size_t target_internal_bytes = 64u * 1024u * 1024u;
    size_t repetitions = target_internal_bytes /
        std::max<size_t>(bytes * passes, 1);
    repetitions = std::max<size_t>(3, std::min<size_t>(repetitions, 100000));
    kernel.function(output, input, bytes, tables, multiplier);

    std::vector<double> nanoseconds;
    std::vector<double> cycles;
    for (unsigned sample = 0; sample < kSamples; ++sample)
    {
        const uint64_t cycle_start = ReadTSCStart();
        const auto time_start = std::chrono::steady_clock::now();
        for (size_t repeat = 0; repeat < repetitions; ++repeat)
            kernel.function(output, input, bytes, tables, multiplier);
        const auto time_stop = std::chrono::steady_clock::now();
        const uint64_t cycle_stop = ReadTSCStop();
        nanoseconds.push_back(std::chrono::duration<double, std::nano>(
            time_stop - time_start).count());
        cycles.push_back(static_cast<double>(cycle_stop - cycle_start));
        gSink += output[(sample * 131u) % bytes];
    }
    const double median_ns = Median(nanoseconds);
    std::vector<double> deviations;
    for (double value : nanoseconds)
        deviations.push_back(std::abs(value - median_ns));
    return { median_ns, Median(deviations), Median(cycles), repetitions };
}

constexpr std::array<uint8_t, 12> kChainCoefficients = {
    1, 0x53, 0xca, 0x02, 0x9d, 0x47, 0xe1, 0x7b, 1, 0x35, 0xfe, 0x11
};

void RunChainOnce(
    const NamedKernel& kernel,
    const Tables& tables,
    uint8_t* a,
    uint8_t* b,
    const uint8_t* side,
    size_t bytes)
{
    uint8_t* source = a;
    uint8_t* destination = b;
    for (uint8_t multiplier : kChainCoefficients)
    {
        kernel.function(destination, source, bytes, tables, multiplier);
        for (size_t i = 0; i < bytes; ++i)
            destination[i] ^= side[i];
        std::swap(source, destination);
    }
    if (source != a)
        std::memcpy(a, source, bytes);
}

Measurement MeasureChain(
    const NamedKernel& kernel,
    const Tables& tables,
    uint8_t* a,
    uint8_t* b,
    const uint8_t* initial,
    const uint8_t* side,
    size_t bytes)
{
    const size_t passes = kChainCoefficients.size();
    const size_t target_internal_bytes = 64u * 1024u * 1024u;
    size_t repetitions = target_internal_bytes /
        std::max<size_t>(bytes * passes, 1);
    repetitions = std::max<size_t>(3, std::min<size_t>(repetitions, 20000));
    std::vector<double> nanoseconds;
    std::vector<double> cycles;
    for (unsigned sample = 0; sample < kSamples; ++sample)
    {
        std::memcpy(a, initial, bytes);
        const uint64_t cycle_start = ReadTSCStart();
        const auto time_start = std::chrono::steady_clock::now();
        for (size_t repeat = 0; repeat < repetitions; ++repeat)
            RunChainOnce(kernel, tables, a, b, side, bytes);
        const auto time_stop = std::chrono::steady_clock::now();
        const uint64_t cycle_stop = ReadTSCStop();
        nanoseconds.push_back(std::chrono::duration<double, std::nano>(
            time_stop - time_start).count());
        cycles.push_back(static_cast<double>(cycle_stop - cycle_start));
        gSink += a[(sample * 131u) % bytes];
    }
    const double median_ns = Median(nanoseconds);
    std::vector<double> deviations;
    for (double value : nanoseconds)
        deviations.push_back(std::abs(value - median_ns));
    return { median_ns, Median(deviations), Median(cycles), repetitions };
}

void PrintMeasurement(
    const char* kind,
    size_t bytes,
    int multiplier,
    unsigned passes,
    const NamedKernel& kernel,
    const Measurement& result)
{
    const double input_bytes = static_cast<double>(bytes) * result.repetitions;
    const double seconds = result.median_ns * 1e-9;
    const double gib = input_bytes / (1024.0 * 1024.0 * 1024.0);
    std::cout << kind << ',' << bytes << ',' << multiplier << ',' << passes << ','
              << kernel.name << ',' << result.repetitions << ',' << kSamples << ','
              << std::fixed << std::setprecision(3)
              << result.median_ns << ',' << result.mad_ns << ','
              << result.median_cycles << ','
              << result.median_cycles / input_bytes << ','
              << gib / seconds << ',' << gib * passes / seconds << '\n';
}

void RunBenchmarks(const Tables& tables)
{
    const std::vector<size_t> sizes = {
        64, 256, 1024, 16 * 1024, 64 * 1024, 1024 * 1024, 16 * 1024 * 1024
    };
    const size_t maximum = sizes.back();
    uint8_t* input = static_cast<uint8_t*>(AlignedAllocate(maximum));
    uint8_t* output = static_cast<uint8_t*>(AlignedAllocate(maximum));
    uint8_t* side = static_cast<uint8_t*>(AlignedAllocate(maximum));
    uint8_t* chain_a = static_cast<uint8_t*>(AlignedAllocate(maximum));
    uint8_t* chain_b = static_cast<uint8_t*>(AlignedAllocate(maximum));
    for (size_t i = 0; i < maximum; ++i)
    {
        input[i] = static_cast<uint8_t>(i * 73u + (i >> 8) * 19u + 11u);
        side[i] = static_cast<uint8_t>(i * 29u + 0xa5u);
    }

    std::cout << "kind,bytes,multiplier,passes,kernel,repetitions,samples,"
                 "median_ns,mad_ns,median_cycles,cycles_per_input_byte,"
                 "effective_gib_s,internal_gib_s\n";
    const std::vector<NamedKernel> kernels = Kernels();
    for (size_t bytes : sizes)
        for (const NamedKernel& kernel : kernels)
            PrintMeasurement("mul", bytes, 0x5d, 1, kernel,
                MeasureMultiply(kernel, tables, output, input, bytes, 0x5d, 1));

    for (uint8_t multiplier : { uint8_t(0), uint8_t(1) })
        for (const NamedKernel& kernel : kernels)
            PrintMeasurement("special", 64 * 1024, multiplier, 1, kernel,
                MeasureMultiply(kernel, tables, output, input, 64 * 1024,
                    multiplier, 1));

    for (size_t bytes : {
            size_t(256), size_t(1024), size_t(4 * 1024),
            size_t(16 * 1024), size_t(64 * 1024),
            size_t(256 * 1024), size_t(1024 * 1024) })
    {
        std::vector<uint8_t> reference(bytes);
        std::vector<uint8_t> temporary(bytes);
        std::memcpy(reference.data(), input, bytes);
        RunChainOnce(kernels[0], tables, reference.data(), temporary.data(), side, bytes);
        for (const NamedKernel& kernel : kernels)
        {
            std::vector<uint8_t> actual(bytes);
            std::vector<uint8_t> work(bytes);
            std::memcpy(actual.data(), input, bytes);
            RunChainOnce(kernel, tables, actual.data(), work.data(), side, bytes);
            Require(actual == reference,
                std::string(kernel.name) + " codec-like chain mismatch");

            PrintMeasurement("chain", bytes, -1,
                static_cast<unsigned>(kChainCoefficients.size()), kernel,
                MeasureChain(kernel, tables, chain_a, chain_b, input, side, bytes));
        }
    }
    std::free(input);
    std::free(output);
    std::free(side);
    std::free(chain_a);
    std::free(chain_b);
}

void PrintMetadata(const Tables& tables)
{
#if defined(__linux__)
    std::cerr << "cpu=" << sched_getcpu() << '\n';
#endif
    std::cerr << "compiler=" << __VERSION__ << '\n';
    std::cerr << "features gfni=" << __builtin_cpu_supports("gfni")
              << " avx2=" << __builtin_cpu_supports("avx2")
              << " avx512f=" << __builtin_cpu_supports("avx512f")
              << " avx512bw=" << __builtin_cpu_supports("avx512bw") << '\n';
    std::cerr << "identity_matrix=0x" << std::hex
              << tables.affine_matrix[1] << std::dec << '\n';
    std::cerr << "aes_root_of_0x11d=0x" << std::hex
              << static_cast<unsigned>(tables.aes_root)
              << " to_aes_matrix=0x" << tables.to_aes_matrix
              << " from_aes_matrix=0x" << tables.from_aes_matrix
              << std::dec << '\n';
    std::cerr << "matrix_hash_fnv1a64=0x" << std::hex
              << FNV1a64(tables.affine_matrix.data(),
                    tables.affine_matrix.size() * sizeof(uint64_t)) << std::dec << '\n';
    std::cerr << "multiply_hash_fnv1a64=0x" << std::hex
              << FNV1a64(tables.multiply.data(),
                    tables.multiply.size() * sizeof(tables.multiply[0]))
              << std::dec << '\n';
    for (unsigned output_bit = 0; output_bit < 8; ++output_bit)
    {
        std::cerr << "operand_bits_for_output_" << output_bit << '=';
        for (unsigned input_bit = 0; input_bit < 8; ++input_bit)
            std::cerr << (input_bit ? "," : "")
                      << tables.operand_bit[output_bit][input_bit];
        std::cerr << '\n';
    }
}

} // namespace

int main(int argc, char** argv)
{
    try
    {
        Require(__builtin_cpu_supports("gfni") &&
                __builtin_cpu_supports("avx2") &&
                __builtin_cpu_supports("avx512f") &&
                __builtin_cpu_supports("avx512bw"),
            "this experiment requires GFNI, AVX2, AVX-512F, and AVX-512BW");
        Tables tables;
        BuildScalarTables(tables);
        DiscoverMatrixBitOrder(tables);
        BuildAffineMatrices(tables);
        BuildAESIsomorphism(tables);
        ValidateAll(tables);
        PrintMetadata(tables);
        std::cerr << "exhaustive_cases=" << (256u * 256u)
                  << " affine_forms=3 mul_isomorphism_forms=3 nibble_forms=2"
                     " status=passed\n";
        if (argc == 2 && std::string(argv[1]) == "--validate")
            return 0;
        Require(argc == 1 || (argc == 2 && std::string(argv[1]) == "--benchmark"),
            "usage: gfni_affine [--validate|--benchmark]");
        RunBenchmarks(tables);
        std::cerr << "sink=" << gSink << '\n';
        return 0;
    }
    catch (const std::exception& error)
    {
        std::cerr << "gfni_affine failed: " << error.what() << '\n';
        return 1;
    }
}
