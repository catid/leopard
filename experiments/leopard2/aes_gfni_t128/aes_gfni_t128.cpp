/*
    Experimental whole-transform AES-representation GFNI screen.

    This is deliberately not a production backend.  It tests the exact
    LEGACY_HIGH GF8 K=65/R=65/B=64 packed encoder shape against Leopard's
    current transform core, including the coordinate-axis gather/scatter and
    the Cantor<->AES conversions.  Production integration is considered only
    if this complete kernel clears the project's whole-call promotion gate.
*/

#include "Leopard2Backend.h"
#include "LeopardFF8.h"
#include "leopard2.h"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <stdexcept>
#include <vector>

#include <immintrin.h>

namespace {

static const unsigned kOriginalCount = 65;
static const unsigned kRecoveryCount = 65;
static const unsigned kSide = 128;
static const unsigned kShardBytes = 64;
static const uint8_t kZeroSkewLog = 255;
static const uint64_t kCantorToAESMatrix = UINT64_C(0xefd0aaca822e5cae);
static const uint64_t kAESToCantorMatrix = UINT64_C(0xffb08466dc727ea0);

class AlignedBuffer
{
public:
    explicit AlignedBuffer(size_t bytes, size_t extra = 0)
        : data_(NULL), bytes_(bytes + extra)
    {
        if (posix_memalign(&data_, 64, bytes_) != 0 || !data_)
            throw std::bad_alloc();
    }

    ~AlignedBuffer()
    {
        std::free(data_);
    }

    uint8_t* bytes()
    {
        return static_cast<uint8_t*>(data_);
    }

private:
    AlignedBuffer(const AlignedBuffer&);
    AlignedBuffer& operator=(const AlignedBuffer&);

    void* data_;
    size_t bytes_;
};

struct Field
{
    std::array<uint8_t, 256> log;
    std::array<uint8_t, 256> exp;
    std::array<uint8_t, 256> cantor_to_aes;
    std::array<uint8_t, 256> aes_to_cantor;
    std::array<uint64_t, 256> affine_matrix;
    std::array<uint8_t, 256> skew_storage;

    static uint8_t AddLog(uint8_t a, uint8_t b)
    {
        const unsigned sum = static_cast<unsigned>(a) + b;
        return static_cast<uint8_t>(sum + (sum >> 8));
    }

    uint8_t MultiplyLog(uint8_t value, uint8_t multiplier_log) const
    {
        if (value == 0)
            return 0;
        return exp[AddLog(log[value], multiplier_log)];
    }

    uint8_t Multiply(uint8_t a, uint8_t b) const
    {
        if (a == 0 || b == 0)
            return 0;
        return exp[AddLog(log[a], log[b])];
    }
};

static uint8_t AESMultiply(uint8_t a, uint8_t b)
{
    unsigned x = a;
    unsigned y = b;
    unsigned result = 0;
    while (y != 0)
    {
        if ((y & 1U) != 0)
            result ^= x;
        y >>= 1;
        x <<= 1;
        if ((x & 0x100U) != 0)
            x ^= 0x11bU;
    }
    return static_cast<uint8_t>(result);
}

static uint8_t AESPower(uint8_t value, unsigned exponent)
{
    uint8_t result = 1;
    while (exponent != 0)
    {
        if ((exponent & 1U) != 0)
            result = AESMultiply(result, value);
        value = AESMultiply(value, value);
        exponent >>= 1;
    }
    return result;
}

static uint8_t CantorToPolynomial(uint8_t value)
{
    static const uint8_t kCantorBasis[8] = {
        1, 214, 152, 146, 86, 200, 88, 230
    };
    uint8_t result = 0;
    for (unsigned bit = 0; bit < 8; ++bit)
        if ((value & (1U << bit)) != 0)
            result ^= kCantorBasis[bit];
    return result;
}

static uint8_t ApplyAffineScalar(uint8_t value, uint64_t matrix)
{
    uint8_t result = 0;
    for (unsigned output_bit = 0; output_bit < 8; ++output_bit)
    {
        const uint8_t row = static_cast<uint8_t>(
            matrix >> (8U * (7U - output_bit)));
        const unsigned parity = static_cast<unsigned>(
            __builtin_parity(static_cast<unsigned>(row & value)));
        result |= static_cast<uint8_t>(parity << output_bit);
    }
    return result;
}

static Field BuildField()
{
    Field field = {};

    unsigned state = 1;
    for (unsigned i = 0; i < 255; ++i)
    {
        field.exp[state] = static_cast<uint8_t>(i);
        state <<= 1;
        if (state >= 256)
            state ^= 0x11dU;
    }
    field.exp[0] = 255;

    static const uint8_t kCantorBasis[8] = {
        1, 214, 152, 146, 86, 200, 88, 230
    };
    field.log[0] = 0;
    for (unsigned bit = 0; bit < 8; ++bit)
    {
        const unsigned width = 1U << bit;
        for (unsigned i = 0; i < width; ++i)
            field.log[i + width] = field.log[i] ^ kCantorBasis[bit];
    }
    for (unsigned i = 0; i < 256; ++i)
        field.log[i] = field.exp[field.log[i]];
    for (unsigned i = 0; i < 256; ++i)
        field.exp[field.log[i]] = static_cast<uint8_t>(i);
    field.exp[255] = field.exp[0];

    const uint8_t aes_root = 0x03;
    for (unsigned cantor = 0; cantor < 256; ++cantor)
    {
        const uint8_t polynomial = CantorToPolynomial(
            static_cast<uint8_t>(cantor));
        uint8_t mapped = 0;
        for (unsigned bit = 0; bit < 8; ++bit)
            if ((polynomial & (1U << bit)) != 0)
                mapped ^= AESPower(aes_root, bit);
        field.cantor_to_aes[cantor] = mapped;
        field.aes_to_cantor[mapped] = static_cast<uint8_t>(cantor);
    }

    uint8_t temp[7];
    for (unsigned i = 1; i < 8; ++i)
        temp[i - 1] = static_cast<uint8_t>(1U << i);
    uint8_t* const skew = field.skew_storage.data() + 1;
    for (unsigned m = 0; m < 7; ++m)
    {
        const unsigned step = 1U << (m + 1U);
        skew[(1U << m) - 1U] = 0;
        for (unsigned i = m; i < 7; ++i)
        {
            const unsigned width = 1U << (i + 1U);
            for (unsigned j = (1U << m) - 1U; j < width; j += step)
                skew[j + width] = skew[j] ^ temp[i];
        }
        temp[m] = static_cast<uint8_t>(255U - field.log[
            field.MultiplyLog(temp[m], field.log[temp[m] ^ 1U])]);
        for (unsigned i = m + 1U; i < 7; ++i)
        {
            const uint8_t sum = Field::AddLog(
                field.log[temp[i] ^ 1U], temp[m]);
            temp[i] = field.MultiplyLog(temp[i], sum);
        }
    }
    for (unsigned i = 0; i < 255; ++i)
        skew[i] = field.log[skew[i]];

    for (unsigned value = 0; value < 256; ++value)
    {
        const uint8_t byte = static_cast<uint8_t>(value);
        if (ApplyAffineScalar(byte, kCantorToAESMatrix) !=
                field.cantor_to_aes[value] ||
            ApplyAffineScalar(field.cantor_to_aes[value],
                kAESToCantorMatrix) != byte)
        {
            throw std::runtime_error("Cantor/AES affine matrix mismatch");
        }
    }
    for (unsigned multiplier = 0; multiplier < 256; ++multiplier)
    {
        uint64_t matrix = 0;
        for (unsigned output_bit = 0; output_bit < 8; ++output_bit)
        {
            uint8_t row = 0;
            for (unsigned input_bit = 0; input_bit < 8; ++input_bit)
            {
                const uint8_t product = field.Multiply(
                    static_cast<uint8_t>(multiplier),
                    static_cast<uint8_t>(1U << input_bit));
                if ((product & (1U << output_bit)) != 0)
                    row |= static_cast<uint8_t>(1U << input_bit);
            }
            matrix |= static_cast<uint64_t>(row) <<
                (8U * (7U - output_bit));
        }
        field.affine_matrix[multiplier] = matrix;
        for (unsigned input = 0; input < 256; ++input)
        {
            const uint8_t expected = field.Multiply(
                static_cast<uint8_t>(multiplier),
                static_cast<uint8_t>(input));
            if (ApplyAffineScalar(static_cast<uint8_t>(input), matrix) !=
                expected)
            {
                throw std::runtime_error(
                    "direct affine multiplier matrix mismatch");
            }
        }
    }
    for (unsigned a = 0; a < 256; ++a)
        for (unsigned b = 0; b < 256; ++b)
            if (AESMultiply(field.cantor_to_aes[a],
                    field.cantor_to_aes[b]) !=
                field.cantor_to_aes[field.Multiply(
                    static_cast<uint8_t>(a), static_cast<uint8_t>(b))])
            {
                throw std::runtime_error("field isomorphism mismatch");
            }

    return field;
}

struct alignas(64) GammaTables
{
    alignas(64) uint8_t inverse[7][128];
    alignas(64) uint8_t forward[7][128];
};

static GammaTables BuildGammaTables(const Field& field)
{
    GammaTables tables = {};
    for (unsigned level = 0; level < 7; ++level)
    {
        const unsigned distance = 1U << level;
        const unsigned block_mask = 2U * distance - 1U;
        for (unsigned coordinate = 0; coordinate < 128; ++coordinate)
        {
            const unsigned block = coordinate & ~block_mask;
            const unsigned skew_index = block + distance;
            const uint8_t inverse_log =
                field.skew_storage[128U + skew_index];
            const uint8_t forward_log = field.skew_storage[skew_index];
            const uint8_t inverse_cantor = inverse_log == kZeroSkewLog
                ? 0 : field.exp[inverse_log];
            const uint8_t forward_cantor = forward_log == kZeroSkewLog
                ? 0 : field.exp[forward_log];
            tables.inverse[level][coordinate] =
                field.cantor_to_aes[inverse_cantor];
            tables.forward[level][coordinate] =
                field.cantor_to_aes[forward_cantor];
        }
    }

    for (unsigned coordinate = 0; coordinate < 128; ++coordinate)
    {
        if (tables.forward[6][coordinate] != 0)
            throw std::runtime_error("forward level-six gamma is not zero");
    }
    for (unsigned coordinate = 0; coordinate < 64; ++coordinate)
    {
        if (tables.forward[5][coordinate] != 0)
            throw std::runtime_error(
                "forward level-five lower gamma is not zero");
    }
    if (tables.inverse[6][0] != 0xbdU)
        throw std::runtime_error("unexpected inverse level-six gamma");
    return tables;
}

static const Field g_field = BuildField();
static const GammaTables g_gamma = BuildGammaTables(g_field);

alignas(64) static const uint8_t kLaneIndices[64] = {
     0,  1,  2,  3,  4,  5,  6,  7,
     8,  9, 10, 11, 12, 13, 14, 15,
    16, 17, 18, 19, 20, 21, 22, 23,
    24, 25, 26, 27, 28, 29, 30, 31,
    32, 33, 34, 35, 36, 37, 38, 39,
    40, 41, 42, 43, 44, 45, 46, 47,
    48, 49, 50, 51, 52, 53, 54, 55,
    56, 57, 58, 59, 60, 61, 62, 63
};

alignas(64) static const int32_t kGatherOffsets[16] = {
      0,   64,  128,  192,  256,  320,  384,  448,
    512,  576,  640,  704,  768,  832,  896,  960
};

#if defined(__GNUC__) || defined(__clang__)
# define AES_TARGET __attribute__((target( \
    "avx2,avx512f,avx512bw,avx512vbmi,gfni")))
# define AES_INLINE inline __attribute__((always_inline, target( \
    "avx2,avx512f,avx512bw,avx512vbmi,gfni")))
#else
# error "This experiment requires GCC or Clang target attributes"
#endif

struct State
{
    __m512i lo;
    __m512i hi;
};

template<unsigned Distance> struct StageMask;
template<> struct StageMask<1> {
    static const uint64_t value = UINT64_C(0xaaaaaaaaaaaaaaaa); };
template<> struct StageMask<2> {
    static const uint64_t value = UINT64_C(0xcccccccccccccccc); };
template<> struct StageMask<4> {
    static const uint64_t value = UINT64_C(0xf0f0f0f0f0f0f0f0); };
template<> struct StageMask<8> {
    static const uint64_t value = UINT64_C(0xff00ff00ff00ff00); };
template<> struct StageMask<16> {
    static const uint64_t value = UINT64_C(0xffff0000ffff0000); };
template<> struct StageMask<32> {
    static const uint64_t value = UINT64_C(0xffffffff00000000); };

template<unsigned Distance>
AES_INLINE __m512i PartnerIndex()
{
    const __m512i lanes = _mm512_load_si512(
        reinterpret_cast<const __m512i*>(kLaneIndices));
    return _mm512_xor_si512(
        lanes, _mm512_set1_epi8(static_cast<char>(Distance)));
}

template<unsigned Distance>
AES_INLINE __m512i InverseStage(
    __m512i value, __m512i partner_index, __m512i gamma)
{
    const __m512i partner =
        _mm512_permutexvar_epi8(partner_index, value);
    const __m512i pair_sum = _mm512_xor_si512(value, partner);
    const __m512i low = _mm512_xor_si512(
        value, _mm512_gf2p8mul_epi8(pair_sum, gamma));
    return _mm512_mask_blend_epi8(
        static_cast<__mmask64>(StageMask<Distance>::value),
        low, pair_sum);
}

template<unsigned Distance>
AES_INLINE __m512i ForwardStage(
    __m512i value, __m512i partner_index, __m512i gamma)
{
    const __m512i partner =
        _mm512_permutexvar_epi8(partner_index, value);
    const __mmask64 high_mask =
        static_cast<__mmask64>(StageMask<Distance>::value);
    const __m512i upper = _mm512_mask_blend_epi8(
        high_mask, partner, value);
    const __m512i low = _mm512_xor_si512(
        value, _mm512_gf2p8mul_epi8(upper, gamma));
    return _mm512_mask_blend_epi8(
        high_mask, low, _mm512_xor_si512(low, partner));
}

template<unsigned Distance>
AES_INLINE __m512i ForwardStageZero(
    __m512i value, __m512i partner_index)
{
    const __m512i partner =
        _mm512_permutexvar_epi8(partner_index, value);
    const __m512i high = _mm512_xor_si512(value, partner);
    return _mm512_mask_blend_epi8(
        static_cast<__mmask64>(StageMask<Distance>::value),
        value, high);
}

template<unsigned Distance>
AES_INLINE void ApplyInverseStage(
    State& a, State& b, State& c, State& d, const uint8_t* gamma)
{
    const __m512i index = PartnerIndex<Distance>();
    const __m512i gamma_lo = _mm512_load_si512(
        reinterpret_cast<const __m512i*>(gamma));
    const __m512i gamma_hi = _mm512_load_si512(
        reinterpret_cast<const __m512i*>(gamma + 64));
    a.lo = InverseStage<Distance>(a.lo, index, gamma_lo);
    a.hi = InverseStage<Distance>(a.hi, index, gamma_hi);
    b.lo = InverseStage<Distance>(b.lo, index, gamma_lo);
    b.hi = InverseStage<Distance>(b.hi, index, gamma_hi);
    c.lo = InverseStage<Distance>(c.lo, index, gamma_lo);
    c.hi = InverseStage<Distance>(c.hi, index, gamma_hi);
    d.lo = InverseStage<Distance>(d.lo, index, gamma_lo);
    d.hi = InverseStage<Distance>(d.hi, index, gamma_hi);
}

template<unsigned Distance>
AES_INLINE void ApplyForwardStage(
    State& a, State& b, State& c, State& d, const uint8_t* gamma)
{
    const __m512i index = PartnerIndex<Distance>();
    const __m512i gamma_lo = _mm512_load_si512(
        reinterpret_cast<const __m512i*>(gamma));
    const __m512i gamma_hi = _mm512_load_si512(
        reinterpret_cast<const __m512i*>(gamma + 64));
    a.lo = ForwardStage<Distance>(a.lo, index, gamma_lo);
    a.hi = ForwardStage<Distance>(a.hi, index, gamma_hi);
    b.lo = ForwardStage<Distance>(b.lo, index, gamma_lo);
    b.hi = ForwardStage<Distance>(b.hi, index, gamma_hi);
    c.lo = ForwardStage<Distance>(c.lo, index, gamma_lo);
    c.hi = ForwardStage<Distance>(c.hi, index, gamma_hi);
    d.lo = ForwardStage<Distance>(d.lo, index, gamma_lo);
    d.hi = ForwardStage<Distance>(d.hi, index, gamma_hi);
}

AES_INLINE void ApplyForwardStage32(
    State& a, State& b, State& c, State& d, const uint8_t* gamma)
{
    const __m512i index = PartnerIndex<32>();
    const __m512i gamma_hi = _mm512_load_si512(
        reinterpret_cast<const __m512i*>(gamma + 64));
    a.lo = ForwardStageZero<32>(a.lo, index);
    a.hi = ForwardStage<32>(a.hi, index, gamma_hi);
    b.lo = ForwardStageZero<32>(b.lo, index);
    b.hi = ForwardStage<32>(b.hi, index, gamma_hi);
    c.lo = ForwardStageZero<32>(c.lo, index);
    c.hi = ForwardStage<32>(c.hi, index, gamma_hi);
    d.lo = ForwardStageZero<32>(d.lo, index);
    d.hi = ForwardStage<32>(d.hi, index, gamma_hi);
}

AES_INLINE void TransformFourStates(
    State& a, State& b, State& c, State& d)
{
    ApplyInverseStage<1>(a, b, c, d, g_gamma.inverse[0]);
    ApplyInverseStage<2>(a, b, c, d, g_gamma.inverse[1]);
    ApplyInverseStage<4>(a, b, c, d, g_gamma.inverse[2]);
    ApplyInverseStage<8>(a, b, c, d, g_gamma.inverse[3]);
    ApplyInverseStage<16>(a, b, c, d, g_gamma.inverse[4]);
    ApplyInverseStage<32>(a, b, c, d, g_gamma.inverse[5]);

    const __m512i inverse_gamma = _mm512_load_si512(
        reinterpret_cast<const __m512i*>(g_gamma.inverse[6]));
#define APPLY_INVERSE_CROSS(state) do { \
        (state).hi = _mm512_xor_si512((state).hi, (state).lo); \
        (state).lo = _mm512_xor_si512((state).lo, \
            _mm512_gf2p8mul_epi8((state).hi, inverse_gamma)); \
    } while (0)
    APPLY_INVERSE_CROSS(a);
    APPLY_INVERSE_CROSS(b);
    APPLY_INVERSE_CROSS(c);
    APPLY_INVERSE_CROSS(d);
#undef APPLY_INVERSE_CROSS

    // The zero-shift forward level-six multiplier is zero.
    a.hi = _mm512_xor_si512(a.hi, a.lo);
    b.hi = _mm512_xor_si512(b.hi, b.lo);
    c.hi = _mm512_xor_si512(c.hi, c.lo);
    d.hi = _mm512_xor_si512(d.hi, d.lo);

    ApplyForwardStage32(a, b, c, d, g_gamma.forward[5]);
    ApplyForwardStage<16>(a, b, c, d, g_gamma.forward[4]);
    ApplyForwardStage<8>(a, b, c, d, g_gamma.forward[3]);
    ApplyForwardStage<4>(a, b, c, d, g_gamma.forward[2]);
    ApplyForwardStage<2>(a, b, c, d, g_gamma.forward[1]);
    ApplyForwardStage<1>(a, b, c, d, g_gamma.forward[0]);
}

template<unsigned ByteIndex>
AES_INLINE State MakeState(
    __m512i rows0, __m512i rows1, __m512i rows2, __m512i rows3,
    uint32_t tail)
{
    const __m128i part0 = _mm512_cvtepi32_epi8(
        _mm512_srli_epi32(rows0, ByteIndex * 8U));
    const __m128i part1 = _mm512_cvtepi32_epi8(
        _mm512_srli_epi32(rows1, ByteIndex * 8U));
    const __m128i part2 = _mm512_cvtepi32_epi8(
        _mm512_srli_epi32(rows2, ByteIndex * 8U));
    const __m128i part3 = _mm512_cvtepi32_epi8(
        _mm512_srli_epi32(rows3, ByteIndex * 8U));
    State state;
    state.lo = _mm512_castsi128_si512(part0);
    state.lo = _mm512_inserti32x4(state.lo, part1, 1);
    state.lo = _mm512_inserti32x4(state.lo, part2, 2);
    state.lo = _mm512_inserti32x4(state.lo, part3, 3);
    state.hi = _mm512_maskz_set1_epi8(
        static_cast<__mmask64>(1),
        static_cast<char>(tail >> (ByteIndex * 8U)));
    return state;
}

AES_INLINE void GatherStates(
    const uint8_t* input, unsigned column,
    State& a, State& b, State& c, State& d)
{
    const __m512i offsets = _mm512_load_si512(
        reinterpret_cast<const __m512i*>(kGatherOffsets));
    const __m512i rows0 = _mm512_i32gather_epi32(
        offsets, input + column, 1);
    const __m512i rows1 = _mm512_i32gather_epi32(
        offsets, input + column + 16U * kShardBytes, 1);
    const __m512i rows2 = _mm512_i32gather_epi32(
        offsets, input + column + 32U * kShardBytes, 1);
    const __m512i rows3 = _mm512_i32gather_epi32(
        offsets, input + column + 48U * kShardBytes, 1);
    uint32_t tail = 0;
    std::memcpy(&tail,
        input + 64U * kShardBytes + column, sizeof(tail));
    a = MakeState<0>(rows0, rows1, rows2, rows3, tail);
    b = MakeState<1>(rows0, rows1, rows2, rows3, tail);
    c = MakeState<2>(rows0, rows1, rows2, rows3, tail);
    d = MakeState<3>(rows0, rows1, rows2, rows3, tail);
}

template<unsigned Block>
AES_INLINE __m512i WidenBlock(__m512i value)
{
    return _mm512_cvtepu8_epi32(
        _mm512_extracti32x4_epi32(value, Block));
}

template<unsigned Block>
AES_INLINE void ScatterBlock(
    uint8_t* output, unsigned column,
    const State& a, const State& b, const State& c, const State& d)
{
    const __m512i a32 = WidenBlock<Block>(a.lo);
    const __m512i b32 = _mm512_slli_epi32(WidenBlock<Block>(b.lo), 8);
    const __m512i c32 = _mm512_slli_epi32(WidenBlock<Block>(c.lo), 16);
    const __m512i d32 = _mm512_slli_epi32(WidenBlock<Block>(d.lo), 24);
    const __m512i packed = _mm512_or_si512(
        _mm512_or_si512(a32, b32), _mm512_or_si512(c32, d32));
    const __m512i offsets = _mm512_load_si512(
        reinterpret_cast<const __m512i*>(kGatherOffsets));
    _mm512_i32scatter_epi32(
        output + column + Block * 16U * kShardBytes,
        offsets, packed, 1);
}

AES_INLINE uint8_t FirstByte(__m512i value)
{
    return static_cast<uint8_t>(_mm_cvtsi128_si32(
        _mm512_castsi512_si128(value)));
}

AES_INLINE void ScatterStates(
    uint8_t* output, unsigned column,
    const State& a, const State& b, const State& c, const State& d)
{
    ScatterBlock<0>(output, column, a, b, c, d);
    ScatterBlock<1>(output, column, a, b, c, d);
    ScatterBlock<2>(output, column, a, b, c, d);
    ScatterBlock<3>(output, column, a, b, c, d);
    const uint32_t tail =
        static_cast<uint32_t>(FirstByte(a.hi)) |
        (static_cast<uint32_t>(FirstByte(b.hi)) << 8) |
        (static_cast<uint32_t>(FirstByte(c.hi)) << 16) |
        (static_cast<uint32_t>(FirstByte(d.hi)) << 24);
    std::memcpy(output + 64U * kShardBytes + column,
        &tail, sizeof(tail));
}

AES_TARGET __attribute__((noinline, aligned(64)))
static void EncodeAESGFNIT128(
    const uint8_t* input, uint8_t* output)
{
    const __m512i to_aes = _mm512_set1_epi64(
        static_cast<long long>(kCantorToAESMatrix));
    const __m512i from_aes = _mm512_set1_epi64(
        static_cast<long long>(kAESToCantorMatrix));
    for (unsigned column = 0; column < kShardBytes; column += 4)
    {
        State a;
        State b;
        State c;
        State d;
        GatherStates(input, column, a, b, c, d);
#define CONVERT_STATE(state, matrix) do { \
        (state).lo = _mm512_gf2p8affine_epi64_epi8( \
            (state).lo, (matrix), 0); \
        (state).hi = _mm512_gf2p8affine_epi64_epi8( \
            (state).hi, (matrix), 0); \
    } while (0)
        CONVERT_STATE(a, to_aes);
        CONVERT_STATE(b, to_aes);
        CONVERT_STATE(c, to_aes);
        CONVERT_STATE(d, to_aes);
        TransformFourStates(a, b, c, d);
        CONVERT_STATE(a, from_aes);
        CONVERT_STATE(b, from_aes);
        CONVERT_STATE(c, from_aes);
        CONVERT_STATE(d, from_aes);
#undef CONVERT_STATE
        ScatterStates(output, column, a, b, c, d);
    }
}

/*
    The coordinate-axis kernel above is the closest translation of the
    supplied report: one transform occupies two ZMM registers.  This second
    screen keeps Leopard's native row-major geometry instead.  Each ZMM is one
    complete 64-byte shard row, so it needs no transpose.  Two product policies
    allow an apples-to-apples comparison between keeping the whole transform
    in AES representation (MULB) and Leopard's direct Cantor-basis affine
    multiplier widened from YMM to ZMM.
*/
AES_INLINE uint8_t CantorGamma(uint8_t multiplier_log)
{
    return multiplier_log == kZeroSkewLog
        ? 0 : g_field.exp[multiplier_log];
}

struct AESRowProduct
{
    typedef uint8_t Constant;

    static AES_INLINE Constant FromLog(uint8_t multiplier_log)
    {
        return g_field.cantor_to_aes[CantorGamma(multiplier_log)];
    }

    static AES_INLINE __m512i Multiply(__m512i value, Constant multiplier)
    {
        return _mm512_gf2p8mul_epi8(
            value, _mm512_set1_epi8(static_cast<char>(multiplier)));
    }

    static AES_INLINE __m512i Enter(__m512i value)
    {
        return _mm512_gf2p8affine_epi64_epi8(
            value,
            _mm512_set1_epi64(
                static_cast<long long>(kCantorToAESMatrix)),
            0);
    }

    static AES_INLINE __m512i Leave(__m512i value)
    {
        return _mm512_gf2p8affine_epi64_epi8(
            value,
            _mm512_set1_epi64(
                static_cast<long long>(kAESToCantorMatrix)),
            0);
    }
};

struct AffineRowProduct
{
    typedef uint64_t Constant;

    static AES_INLINE Constant FromLog(uint8_t multiplier_log)
    {
        return g_field.affine_matrix[CantorGamma(multiplier_log)];
    }

    static AES_INLINE __m512i Multiply(__m512i value, Constant multiplier)
    {
        return _mm512_gf2p8affine_epi64_epi8(
            value,
            _mm512_set1_epi64(static_cast<long long>(multiplier)),
            0);
    }

    static AES_INLINE __m512i Enter(__m512i value)
    {
        return value;
    }

    static AES_INLINE __m512i Leave(__m512i value)
    {
        return value;
    }
};

AES_INLINE __m512i LoadRow(const uint8_t* state, unsigned row)
{
    return _mm512_load_si512(reinterpret_cast<const __m512i*>(
        state + row * kShardBytes));
}

AES_INLINE void StoreRow(uint8_t* state, unsigned row, __m512i value)
{
    _mm512_store_si512(reinterpret_cast<__m512i*>(
        state + row * kShardBytes), value);
}

template<class Product>
AES_INLINE void InverseRow4(
    uint8_t* state, unsigned row, unsigned distance,
    typename Product::Constant multiplier01,
    typename Product::Constant multiplier23,
    typename Product::Constant multiplier02)
{
    __m512i x0 = LoadRow(state, row);
    __m512i x1 = LoadRow(state, row + distance);
    __m512i x2 = LoadRow(state, row + distance * 2U);
    __m512i x3 = LoadRow(state, row + distance * 3U);

    x1 = _mm512_xor_si512(x1, x0);
    if (multiplier01 != 0)
        x0 = _mm512_xor_si512(
            x0, Product::Multiply(x1, multiplier01));
    x3 = _mm512_xor_si512(x3, x2);
    if (multiplier23 != 0)
        x2 = _mm512_xor_si512(
            x2, Product::Multiply(x3, multiplier23));
    x2 = _mm512_xor_si512(x2, x0);
    x3 = _mm512_xor_si512(x3, x1);
    if (multiplier02 != 0)
    {
        x0 = _mm512_xor_si512(
            x0, Product::Multiply(x2, multiplier02));
        x1 = _mm512_xor_si512(
            x1, Product::Multiply(x3, multiplier02));
    }

    StoreRow(state, row, x0);
    StoreRow(state, row + distance, x1);
    StoreRow(state, row + distance * 2U, x2);
    StoreRow(state, row + distance * 3U, x3);
}

template<class Product>
AES_INLINE void ForwardRow4(
    uint8_t* state, unsigned row, unsigned distance,
    typename Product::Constant multiplier01,
    typename Product::Constant multiplier23,
    typename Product::Constant multiplier02)
{
    __m512i x0 = LoadRow(state, row);
    __m512i x1 = LoadRow(state, row + distance);
    __m512i x2 = LoadRow(state, row + distance * 2U);
    __m512i x3 = LoadRow(state, row + distance * 3U);

    if (multiplier02 != 0)
    {
        x0 = _mm512_xor_si512(
            x0, Product::Multiply(x2, multiplier02));
        x1 = _mm512_xor_si512(
            x1, Product::Multiply(x3, multiplier02));
    }
    x2 = _mm512_xor_si512(x2, x0);
    x3 = _mm512_xor_si512(x3, x1);
    if (multiplier01 != 0)
        x0 = _mm512_xor_si512(
            x0, Product::Multiply(x1, multiplier01));
    x1 = _mm512_xor_si512(x1, x0);
    if (multiplier23 != 0)
        x2 = _mm512_xor_si512(
            x2, Product::Multiply(x3, multiplier23));
    x3 = _mm512_xor_si512(x3, x2);

    StoreRow(state, row, x0);
    StoreRow(state, row + distance, x1);
    StoreRow(state, row + distance * 2U, x2);
    StoreRow(state, row + distance * 3U, x3);
}

template<class Product>
AES_INLINE void InverseRow2(
    uint8_t* state, unsigned row, unsigned distance,
    typename Product::Constant multiplier)
{
    __m512i x = LoadRow(state, row);
    __m512i y = LoadRow(state, row + distance);
    y = _mm512_xor_si512(y, x);
    if (multiplier != 0)
        x = _mm512_xor_si512(x, Product::Multiply(y, multiplier));
    StoreRow(state, row, x);
    StoreRow(state, row + distance, y);
}

template<class Product>
AES_INLINE void ForwardRow2(
    uint8_t* state, unsigned row, unsigned distance,
    typename Product::Constant multiplier)
{
    __m512i x = LoadRow(state, row);
    __m512i y = LoadRow(state, row + distance);
    if (multiplier != 0)
        x = _mm512_xor_si512(x, Product::Multiply(y, multiplier));
    y = _mm512_xor_si512(y, x);
    StoreRow(state, row, x);
    StoreRow(state, row + distance, y);
}

template<class Product>
AES_INLINE void RowTransform(
    const uint8_t* input, uint8_t* output, uint8_t* state)
{
    for (unsigned row = 0; row < kOriginalCount; ++row)
    {
        const __m512i value = _mm512_loadu_si512(
            input + row * kShardBytes);
        StoreRow(state, row, Product::Enter(value));
    }
    const __m512i zero = _mm512_setzero_si512();
    for (unsigned row = kOriginalCount; row < kSide; ++row)
        StoreRow(state, row, zero);

    for (unsigned distance = 1; distance <= 16; distance *= 4)
    {
        const unsigned width = distance * 4U;
        for (unsigned row = 0; row < kOriginalCount; row += width)
        {
            const uint8_t* const skew =
                g_field.skew_storage.data() + 128U + row;
            const typename Product::Constant multiplier01 =
                Product::FromLog(skew[distance]);
            const typename Product::Constant multiplier23 =
                Product::FromLog(skew[distance * 3U]);
            const typename Product::Constant multiplier02 =
                Product::FromLog(skew[distance * 2U]);
            for (unsigned lane = 0; lane < distance; ++lane)
            {
                InverseRow4<Product>(state, row + lane, distance,
                    multiplier01, multiplier23, multiplier02);
            }
        }
    }
    const typename Product::Constant inverse_final =
        Product::FromLog(g_field.skew_storage[128U + 64U]);
    for (unsigned row = 0; row < 64; ++row)
        InverseRow2<Product>(state, row, 64, inverse_final);

    for (unsigned distance = 32; distance >= 2; distance /= 4)
    {
        const unsigned width = distance * 4U;
        for (unsigned row = 0; row < kRecoveryCount; row += width)
        {
            const uint8_t* const skew =
                g_field.skew_storage.data() + row;
            const typename Product::Constant multiplier01 =
                Product::FromLog(skew[distance]);
            const typename Product::Constant multiplier23 =
                Product::FromLog(skew[distance * 3U]);
            const typename Product::Constant multiplier02 =
                Product::FromLog(skew[distance * 2U]);
            for (unsigned lane = 0; lane < distance; ++lane)
            {
                ForwardRow4<Product>(state, row + lane, distance,
                    multiplier01, multiplier23, multiplier02);
            }
        }
    }
    for (unsigned row = 0; row < kRecoveryCount; row += 2)
    {
        const typename Product::Constant multiplier =
            Product::FromLog(g_field.skew_storage[row + 1U]);
        ForwardRow2<Product>(state, row, 1, multiplier);
    }

    for (unsigned row = 0; row < kRecoveryCount; ++row)
    {
        const __m512i value = Product::Leave(LoadRow(state, row));
        _mm512_storeu_si512(output + row * kShardBytes, value);
    }
}

AES_TARGET __attribute__((noinline, aligned(64)))
static void EncodeRowAESGFNI(
    const uint8_t* input, uint8_t* output, uint8_t* state)
{
    RowTransform<AESRowProduct>(input, output, state);
}

AES_TARGET __attribute__((noinline, aligned(64)))
static void EncodeRowAffineGFNI(
    const uint8_t* input, uint8_t* output, uint8_t* state)
{
    RowTransform<AffineRowProduct>(input, output, state);
}

static bool HostSupportsCandidate()
{
#if defined(__GNUC__) || defined(__clang__)
    __builtin_cpu_init();
    return __builtin_cpu_supports("avx2") &&
        __builtin_cpu_supports("avx512f") &&
        __builtin_cpu_supports("avx512bw") &&
        __builtin_cpu_supports("avx512vbmi") &&
        __builtin_cpu_supports("gfni");
#else
    return false;
#endif
}

static uint64_t SplitMix64(uint64_t& state)
{
    uint64_t z = (state += UINT64_C(0x9e3779b97f4a7c15));
    z = (z ^ (z >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
    z = (z ^ (z >> 27)) * UINT64_C(0x94d049bb133111eb);
    return z ^ (z >> 31);
}

static void FillInput(uint8_t* input, uint64_t seed)
{
    for (unsigned i = 0; i < kOriginalCount * kShardBytes; ++i)
        input[i] = static_cast<uint8_t>(SplitMix64(seed));
}

static void RunReference(
    const leopard::backend::Ops& ops,
    const void* const* data,
    void** work)
{
    leopard::ff8::ReedSolomonEncodeK65R65T128B64Packed(
        ops, data, work);
}

static void RequireEqual(
    const uint8_t* expected, const uint8_t* actual, const char* context)
{
    const size_t bytes = kRecoveryCount * kShardBytes;
    for (size_t i = 0; i < bytes; ++i)
        if (expected[i] != actual[i])
        {
            std::fprintf(stderr,
                "%s mismatch at byte %zu: expected=%u actual=%u\n",
                context, i, static_cast<unsigned>(expected[i]),
                static_cast<unsigned>(actual[i]));
            throw std::runtime_error("candidate parity mismatch");
        }
}

static void VerifyCandidate(
    const leopard::backend::Ops& ops,
    uint8_t* input,
    uint8_t* coordinate_candidate,
    uint8_t* row_aes_candidate,
    uint8_t* row_affine_candidate,
    uint8_t* row_state,
    uint8_t* reference_storage,
    const void* data[65],
    void* work[128])
{
    FillInput(input, UINT64_C(0x6b36357236356236));
    RunReference(ops, data, work);
    EncodeAESGFNIT128(input, coordinate_candidate);
    EncodeRowAESGFNI(input, row_aes_candidate, row_state);
    EncodeRowAffineGFNI(input, row_affine_candidate, row_state);
    RequireEqual(reference_storage, coordinate_candidate,
        "random coordinate-axis input");
    RequireEqual(reference_storage, row_aes_candidate,
        "random row-major AES input");
    RequireEqual(reference_storage, row_affine_candidate,
        "random row-major affine input");

    for (unsigned source = 0; source < kOriginalCount; ++source)
    {
        std::memset(input, 0, kOriginalCount * kShardBytes);
        for (unsigned column = 0; column < kShardBytes; ++column)
            input[source * kShardBytes + column] =
                static_cast<uint8_t>(column * 29U + source * 17U + 1U);
        RunReference(ops, data, work);
        EncodeAESGFNIT128(input, coordinate_candidate);
        EncodeRowAESGFNI(input, row_aes_candidate, row_state);
        EncodeRowAffineGFNI(input, row_affine_candidate, row_state);
        RequireEqual(reference_storage, coordinate_candidate,
            "coordinate-axis systematic basis");
        RequireEqual(reference_storage, row_aes_candidate,
            "row-major AES systematic basis");
        RequireEqual(reference_storage, row_affine_candidate,
            "row-major affine systematic basis");
    }

    for (unsigned iteration = 0; iteration < 128; ++iteration)
    {
        FillInput(input,
            UINT64_C(0xabcdef0123456789) + iteration * UINT64_C(0x101));
        RunReference(ops, data, work);
        EncodeAESGFNIT128(input, coordinate_candidate);
        EncodeRowAESGFNI(input, row_aes_candidate, row_state);
        EncodeRowAffineGFNI(input, row_affine_candidate, row_state);
        RequireEqual(reference_storage, coordinate_candidate,
            "coordinate-axis random differential");
        RequireEqual(reference_storage, row_aes_candidate,
            "row-major AES random differential");
        RequireEqual(reference_storage, row_affine_candidate,
            "row-major affine random differential");
    }
}

static double Median(std::vector<double> values)
{
    std::sort(values.begin(), values.end());
    const size_t middle = values.size() / 2;
    if ((values.size() & 1U) != 0)
        return values[middle];
    return 0.5 * (values[middle - 1] + values[middle]);
}

static double MAD(const std::vector<double>& values, double median)
{
    std::vector<double> deviations;
    deviations.reserve(values.size());
    for (size_t i = 0; i < values.size(); ++i)
        deviations.push_back(std::fabs(values[i] - median));
    return Median(deviations);
}

template<class Function>
static double TimeCalls(Function function, unsigned repeats)
{
    const std::chrono::steady_clock::time_point begin =
        std::chrono::steady_clock::now();
    for (unsigned i = 0; i < repeats; ++i)
        function();
    const std::chrono::steady_clock::time_point end =
        std::chrono::steady_clock::now();
    return std::chrono::duration<double, std::nano>(end - begin).count() /
        repeats;
}

static volatile uint8_t g_sink = 0;

static void Benchmark(
    const leopard::backend::Ops& ops,
    uint8_t* input,
    uint8_t* coordinate_candidate,
    uint8_t* row_aes_candidate,
    uint8_t* row_affine_candidate,
    uint8_t* row_state,
    uint8_t* reference_storage,
    const void* data[65],
    void* work[128])
{
    static const unsigned kWarmup = 2048;
    static const unsigned kSamples = 31;
    static const unsigned kRepeats = 4096;

    FillInput(input, UINT64_C(0x123456789abcdef0));
    for (unsigned i = 0; i < kWarmup; ++i)
    {
        RunReference(ops, data, work);
        EncodeAESGFNIT128(input, coordinate_candidate);
        EncodeRowAESGFNI(input, row_aes_candidate, row_state);
        EncodeRowAffineGFNI(input, row_affine_candidate, row_state);
    }

    std::vector<double> reference_ns;
    std::vector<double> coordinate_ns;
    std::vector<double> row_aes_ns;
    std::vector<double> row_affine_ns;
    std::vector<double> coordinate_ratios;
    std::vector<double> row_aes_ratios;
    std::vector<double> row_affine_ratios;
    reference_ns.reserve(kSamples);
    coordinate_ns.reserve(kSamples);
    row_aes_ns.reserve(kSamples);
    row_affine_ns.reserve(kSamples);
    coordinate_ratios.reserve(kSamples);
    row_aes_ratios.reserve(kSamples);
    row_affine_ratios.reserve(kSamples);
    for (unsigned sample = 0; sample < kSamples; ++sample)
    {
        double reference = 0;
        double coordinate = 0;
        double row_aes = 0;
        double row_affine = 0;
        const auto reference_call = [&]() {
            RunReference(ops, data, work);
            g_sink ^= reference_storage[(sample * 131U) %
                (kRecoveryCount * kShardBytes)];
        };
        const auto coordinate_call = [&]() {
            EncodeAESGFNIT128(input, coordinate_candidate);
            g_sink ^= coordinate_candidate[(sample * 193U) %
                (kRecoveryCount * kShardBytes)];
        };
        const auto row_aes_call = [&]() {
            EncodeRowAESGFNI(input, row_aes_candidate, row_state);
            g_sink ^= row_aes_candidate[(sample * 197U) %
                (kRecoveryCount * kShardBytes)];
        };
        const auto row_affine_call = [&]() {
            EncodeRowAffineGFNI(input, row_affine_candidate, row_state);
            g_sink ^= row_affine_candidate[(sample * 199U) %
                (kRecoveryCount * kShardBytes)];
        };
        switch (sample & 3U)
        {
        case 0:
            reference = TimeCalls(reference_call, kRepeats);
            coordinate = TimeCalls(coordinate_call, kRepeats);
            row_aes = TimeCalls(row_aes_call, kRepeats);
            row_affine = TimeCalls(row_affine_call, kRepeats);
            break;
        case 1:
            row_affine = TimeCalls(row_affine_call, kRepeats);
            row_aes = TimeCalls(row_aes_call, kRepeats);
            coordinate = TimeCalls(coordinate_call, kRepeats);
            reference = TimeCalls(reference_call, kRepeats);
            break;
        case 2:
            row_aes = TimeCalls(row_aes_call, kRepeats);
            reference = TimeCalls(reference_call, kRepeats);
            row_affine = TimeCalls(row_affine_call, kRepeats);
            coordinate = TimeCalls(coordinate_call, kRepeats);
            break;
        default:
            coordinate = TimeCalls(coordinate_call, kRepeats);
            row_affine = TimeCalls(row_affine_call, kRepeats);
            reference = TimeCalls(reference_call, kRepeats);
            row_aes = TimeCalls(row_aes_call, kRepeats);
            break;
        }
        reference_ns.push_back(reference);
        coordinate_ns.push_back(coordinate);
        row_aes_ns.push_back(row_aes);
        row_affine_ns.push_back(row_affine);
        coordinate_ratios.push_back(reference / coordinate);
        row_aes_ratios.push_back(reference / row_aes);
        row_affine_ratios.push_back(reference / row_affine);
    }

    const double reference_median = Median(reference_ns);
    std::printf("reference_core_median_ns=%.3f mad_ns=%.3f\n",
        reference_median, MAD(reference_ns, reference_median));
#define PRINT_RESULT(label, timings, ratios) do { \
        const double timing_median = Median(timings); \
        const double ratio_median = Median(ratios); \
        std::printf(label "_median_ns=%.3f mad_ns=%.3f\n", \
            timing_median, MAD(timings, timing_median)); \
        std::printf(label "_paired_speedup=%.6fx mad=%.6f " \
            "promotion_point_estimate=%s\n", \
            ratio_median, MAD(ratios, ratio_median), \
            ratio_median >= 1.05 ? "pass" : "fail"); \
    } while (0)
    PRINT_RESULT("coordinate_aes_gfni", coordinate_ns, coordinate_ratios);
    PRINT_RESULT("row_aes_gfni", row_aes_ns, row_aes_ratios);
    PRINT_RESULT("row_affine_gfni", row_affine_ns, row_affine_ratios);
#undef PRINT_RESULT
}

} // namespace

int main()
{
    try
    {
        if (!HostSupportsCandidate())
        {
            std::printf(
                "SKIP: AVX2+AVX512F/BW/VBMI+GFNI with OS ZMM state required\n");
            return 0;
        }

        leo2_context_options options = {};
        options.struct_size = sizeof(options);
        options.backend = LEO2_BACKEND_AVX2;
        options.thread_count = 1;
        leo2_context* context = NULL;
        const leo2_result result = leo2_context_create(&options, &context);
        if (result != LEO2_SUCCESS || !context)
            throw std::runtime_error("could not create explicit AVX2 context");

        leopard::backend::QualificationStatus status =
            leopard::backend::QualificationAvailable;
        const leopard::backend::Ops* ops =
            leopard::backend::GetQualifiedOps(LEO2_BACKEND_AVX2, &status);
        if (!ops || ops->kind != LEO2_BACKEND_AVX2)
            throw std::runtime_error("could not qualify AVX2 ops");

        AlignedBuffer input_allocation(
            kOriginalCount * kShardBytes, 8);
        AlignedBuffer coordinate_allocation(
            kRecoveryCount * kShardBytes, 8);
        AlignedBuffer row_aes_allocation(
            kRecoveryCount * kShardBytes, 8);
        AlignedBuffer row_affine_allocation(
            kRecoveryCount * kShardBytes, 8);
        AlignedBuffer row_state_allocation(kSide * kShardBytes);
        AlignedBuffer reference_allocation(kSide * kShardBytes);
        uint8_t* const input = input_allocation.bytes() + 1;
        uint8_t* const coordinate_candidate =
            coordinate_allocation.bytes() + 3;
        uint8_t* const row_aes_candidate =
            row_aes_allocation.bytes() + 5;
        uint8_t* const row_affine_candidate =
            row_affine_allocation.bytes() + 7;
        uint8_t* const row_state = row_state_allocation.bytes();
        uint8_t* const reference_storage = reference_allocation.bytes();
        const void* data[kOriginalCount];
        void* work[kSide];
        for (unsigned row = 0; row < kOriginalCount; ++row)
            data[row] = input + row * kShardBytes;
        for (unsigned row = 0; row < kSide; ++row)
            work[row] = reference_storage + row * kShardBytes;

        VerifyCandidate(*ops, input, coordinate_candidate,
            row_aes_candidate, row_affine_candidate, row_state,
            reference_storage, data, work);
        std::printf("correctness=pass basis_cases=65 random_cases=129\n");
        Benchmark(*ops, input, coordinate_candidate,
            row_aes_candidate, row_affine_candidate, row_state,
            reference_storage, data, work);
        leo2_context_destroy(context);
        return 0;
    }
    catch (const std::exception& error)
    {
        std::fprintf(stderr, "FAIL: %s\n", error.what());
        return 1;
    }
}
