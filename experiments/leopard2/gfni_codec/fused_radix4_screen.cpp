// Radix-four IFFT layer-pair: current nibble split vs GFNI split vs GFNI fused
// with memory-operand matrices.  256-bit data path throughout.
//
// A radix-four inverse group is, in order:
//   butterfly(v0,v1,log01), butterfly(v2,v3,log23),
//   butterfly(v0,v2,log02), butterfly(v1,v3,log02)
// with inverse butterfly(x,y,L): y ^= x; x ^= mul(y,L).
//
// The three variants must agree byte for byte; the harness checks that first.
#include <immintrin.h>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <chrono>
#include <vector>
#include <algorithm>

namespace {

typedef std::chrono::steady_clock Clock;

uint64_t Rng(uint64_t& s)
{
    s ^= s << 13; s ^= s >> 7; s ^= s << 17;
    return s;
}

// ------------------------------------------------------------------ GF8 maps

struct Map8
{
    uint8_t low[16];
    uint8_t high[16];
    uint8_t affine[32];   // matrix replicated four times for a ymm memory operand
};

Map8 MakeMap8(uint64_t& seed)
{
    Map8 map;
    uint8_t column[8];
    for (unsigned i = 0; i < 8; ++i)
        column[i] = static_cast<uint8_t>(Rng(seed));
    for (unsigned n = 0; n < 16; ++n)
    {
        uint8_t lo = 0, hi = 0;
        for (unsigned b = 0; b < 4; ++b)
            if (n & (1u << b)) { lo ^= column[b]; hi ^= column[b + 4]; }
        map.low[n] = lo;
        map.high[n] = hi;
    }
    uint64_t matrix = 0;
    for (unsigned out = 0; out < 8; ++out)
        for (unsigned in = 0; in < 8; ++in)
            if (column[in] & (1u << out))
                matrix |= 1ULL << (8 * (7 - out) + in);
    for (unsigned rep = 0; rep < 4; ++rep)
        memcpy(map.affine + rep * 8, &matrix, 8);
    return map;
}

__attribute__((target("avx2")))
inline __m256i Product8Nibble(__m256i data, __m256i lo, __m256i hi)
{
    const __m256i mask = _mm256_set1_epi8(15);
    return _mm256_xor_si256(
        _mm256_shuffle_epi8(lo, _mm256_and_si256(data, mask)),
        _mm256_shuffle_epi8(hi,
            _mm256_and_si256(_mm256_srli_epi64(data, 4), mask)));
}

__attribute__((target("avx2")))
inline __m256i Bcast(const uint8_t t[16])
{
    return _mm256_broadcastsi128_si256(
        _mm_loadu_si128(reinterpret_cast<const __m128i*>(t)));
}

// (a) today's schedule: four independent radix-two sweeps, nibble tables.
__attribute__((target("avx2")))
void Gf8SplitNibble(uint8_t* v0, uint8_t* v1, uint8_t* v2, uint8_t* v3,
                    const Map8& m01, const Map8& m23, const Map8& m02,
                    size_t bytes)
{
    struct Pair { uint8_t* x; uint8_t* y; const Map8* m; };
    const Pair layers[4] = {
        { v0, v1, &m01 }, { v2, v3, &m23 }, { v0, v2, &m02 }, { v1, v3, &m02 }
    };
    for (unsigned layer = 0; layer < 4; ++layer)
    {
        const __m256i lo = Bcast(layers[layer].m->low);
        const __m256i hi = Bcast(layers[layer].m->high);
        uint8_t* x = layers[layer].x;
        uint8_t* y = layers[layer].y;
        for (size_t i = 0; i < bytes; i += 32)
        {
            __m256i xv = _mm256_loadu_si256((const __m256i*)(x + i));
            __m256i yv = _mm256_loadu_si256((const __m256i*)(y + i));
            yv = _mm256_xor_si256(yv, xv);
            xv = _mm256_xor_si256(xv, Product8Nibble(yv, lo, hi));
            _mm256_storeu_si256((__m256i*)(x + i), xv);
            _mm256_storeu_si256((__m256i*)(y + i), yv);
        }
    }
}

// (b) same schedule, affine multiply in registers.
__attribute__((target("avx2,gfni")))
void Gf8SplitGfni(uint8_t* v0, uint8_t* v1, uint8_t* v2, uint8_t* v3,
                  const Map8& m01, const Map8& m23, const Map8& m02,
                  size_t bytes)
{
    struct Pair { uint8_t* x; uint8_t* y; const Map8* m; };
    const Pair layers[4] = {
        { v0, v1, &m01 }, { v2, v3, &m23 }, { v0, v2, &m02 }, { v1, v3, &m02 }
    };
    for (unsigned layer = 0; layer < 4; ++layer)
    {
        const __m256i mat =
            _mm256_loadu_si256((const __m256i*)layers[layer].m->affine);
        uint8_t* x = layers[layer].x;
        uint8_t* y = layers[layer].y;
        for (size_t i = 0; i < bytes; i += 32)
        {
            __m256i xv = _mm256_loadu_si256((const __m256i*)(x + i));
            __m256i yv = _mm256_loadu_si256((const __m256i*)(y + i));
            yv = _mm256_xor_si256(yv, xv);
            xv = _mm256_xor_si256(xv,
                _mm256_gf2p8affine_epi64_epi8(yv, mat, 0));
            _mm256_storeu_si256((__m256i*)(x + i), xv);
            _mm256_storeu_si256((__m256i*)(y + i), yv);
        }
    }
}

// (c) fused radix-four, matrices as memory operands.
__attribute__((target("avx2,gfni")))
void Gf8FusedGfni(uint8_t* v0, uint8_t* v1, uint8_t* v2, uint8_t* v3,
                  const Map8& m01, const Map8& m23, const Map8& m02,
                  size_t bytes)
{
    const __m256i* p01 = (const __m256i*)m01.affine;
    const __m256i* p23 = (const __m256i*)m23.affine;
    const __m256i* p02 = (const __m256i*)m02.affine;
    for (size_t i = 0; i < bytes; i += 32)
    {
        __m256i a = _mm256_loadu_si256((const __m256i*)(v0 + i));
        __m256i b = _mm256_loadu_si256((const __m256i*)(v1 + i));
        __m256i c = _mm256_loadu_si256((const __m256i*)(v2 + i));
        __m256i d = _mm256_loadu_si256((const __m256i*)(v3 + i));
        b = _mm256_xor_si256(b, a);
        a = _mm256_xor_si256(a, _mm256_gf2p8affine_epi64_epi8(
            b, _mm256_loadu_si256(p01), 0));
        d = _mm256_xor_si256(d, c);
        c = _mm256_xor_si256(c, _mm256_gf2p8affine_epi64_epi8(
            d, _mm256_loadu_si256(p23), 0));
        c = _mm256_xor_si256(c, a);
        a = _mm256_xor_si256(a, _mm256_gf2p8affine_epi64_epi8(
            c, _mm256_loadu_si256(p02), 0));
        d = _mm256_xor_si256(d, b);
        b = _mm256_xor_si256(b, _mm256_gf2p8affine_epi64_epi8(
            d, _mm256_loadu_si256(p02), 0));
        _mm256_storeu_si256((__m256i*)(v0 + i), a);
        _mm256_storeu_si256((__m256i*)(v1 + i), b);
        _mm256_storeu_si256((__m256i*)(v2 + i), c);
        _mm256_storeu_si256((__m256i*)(v3 + i), d);
    }
}


// (a2) production GF8 baseline: fused radix-four with nibble tables in
// registers, mirroring AVX2FF8Butterfly4RangePreparedImpl.
__attribute__((target("avx2")))
void Gf8FusedNibble(uint8_t* v0, uint8_t* v1, uint8_t* v2, uint8_t* v3,
                    const Map8& m01, const Map8& m23, const Map8& m02,
                    size_t bytes)
{
    const __m256i lo01 = Bcast(m01.low), hi01 = Bcast(m01.high);
    const __m256i lo23 = Bcast(m23.low), hi23 = Bcast(m23.high);
    const __m256i lo02 = Bcast(m02.low), hi02 = Bcast(m02.high);
    for (size_t i = 0; i < bytes; i += 32)
    {
        __m256i a = _mm256_loadu_si256((const __m256i*)(v0 + i));
        __m256i b = _mm256_loadu_si256((const __m256i*)(v1 + i));
        __m256i c = _mm256_loadu_si256((const __m256i*)(v2 + i));
        __m256i d = _mm256_loadu_si256((const __m256i*)(v3 + i));
        b = _mm256_xor_si256(b, a);
        a = _mm256_xor_si256(a, Product8Nibble(b, lo01, hi01));
        d = _mm256_xor_si256(d, c);
        c = _mm256_xor_si256(c, Product8Nibble(d, lo23, hi23));
        c = _mm256_xor_si256(c, a);
        a = _mm256_xor_si256(a, Product8Nibble(c, lo02, hi02));
        d = _mm256_xor_si256(d, b);
        b = _mm256_xor_si256(b, Product8Nibble(d, lo02, hi02));
        _mm256_storeu_si256((__m256i*)(v0 + i), a);
        _mm256_storeu_si256((__m256i*)(v1 + i), b);
        _mm256_storeu_si256((__m256i*)(v2 + i), c);
        _mm256_storeu_si256((__m256i*)(v3 + i), d);
    }
}

// ----------------------------------------------------------------- GF16 maps

struct Map16
{
    uint8_t low[4][16];
    uint8_t high[4][16];
    uint8_t affine[4][32];   // A(lo<-lo), A(lo<-hi), A(hi<-lo), A(hi<-hi)
};

Map16 MakeMap16(uint64_t& seed)
{
    Map16 map;
    uint16_t column[16];
    for (unsigned i = 0; i < 16; ++i)
        column[i] = static_cast<uint16_t>(Rng(seed));
    for (unsigned part = 0; part < 4; ++part)
        for (unsigned n = 0; n < 16; ++n)
        {
            uint16_t value = 0;
            for (unsigned b = 0; b < 4; ++b)
                if (n & (1u << b)) value ^= column[part * 4 + b];
            map.low[part][n] = static_cast<uint8_t>(value);
            map.high[part][n] = static_cast<uint8_t>(value >> 8);
        }
    for (unsigned block = 0; block < 4; ++block)
    {
        const unsigned out_half = block >> 1;
        const unsigned in_half = block & 1;
        uint64_t matrix = 0;
        for (unsigned out = 0; out < 8; ++out)
            for (unsigned in = 0; in < 8; ++in)
            {
                const uint16_t product = column[in_half * 8 + in];
                const uint8_t byte = static_cast<uint8_t>(
                    out_half ? (product >> 8) : product);
                if (byte & (1u << out))
                    matrix |= 1ULL << (8 * (7 - out) + in);
            }
        for (unsigned rep = 0; rep < 4; ++rep)
            memcpy(map.affine[block] + rep * 8, &matrix, 8);
    }
    return map;
}

__attribute__((target("avx2")))
inline void Product16Nibble(__m256i dl, __m256i dh,
                            const __m256i lo[4], const __m256i hi[4],
                            __m256i& pl, __m256i& ph)
{
    const __m256i mask = _mm256_set1_epi8(15);
    const __m256i n[4] = {
        _mm256_and_si256(dl, mask),
        _mm256_and_si256(_mm256_srli_epi64(dl, 4), mask),
        _mm256_and_si256(dh, mask),
        _mm256_and_si256(_mm256_srli_epi64(dh, 4), mask)
    };
    pl = _mm256_setzero_si256();
    ph = _mm256_setzero_si256();
    for (unsigned i = 0; i < 4; ++i)
    {
        pl = _mm256_xor_si256(pl, _mm256_shuffle_epi8(lo[i], n[i]));
        ph = _mm256_xor_si256(ph, _mm256_shuffle_epi8(hi[i], n[i]));
    }
}

__attribute__((target("avx2")))
void Gf16SplitNibble(uint8_t* v0, uint8_t* v1, uint8_t* v2, uint8_t* v3,
                     const Map16& m01, const Map16& m23, const Map16& m02,
                     size_t bytes)
{
    struct Pair { uint8_t* x; uint8_t* y; const Map16* m; };
    const Pair layers[4] = {
        { v0, v1, &m01 }, { v2, v3, &m23 }, { v0, v2, &m02 }, { v1, v3, &m02 }
    };
    for (unsigned layer = 0; layer < 4; ++layer)
    {
        __m256i lo[4], hi[4];
        for (unsigned i = 0; i < 4; ++i)
        {
            lo[i] = Bcast(layers[layer].m->low[i]);
            hi[i] = Bcast(layers[layer].m->high[i]);
        }
        uint8_t* x = layers[layer].x;
        uint8_t* y = layers[layer].y;
        for (size_t i = 0; i < bytes; i += 64)
        {
            __m256i xl = _mm256_loadu_si256((const __m256i*)(x + i));
            __m256i xh = _mm256_loadu_si256((const __m256i*)(x + i + 32));
            __m256i yl = _mm256_loadu_si256((const __m256i*)(y + i));
            __m256i yh = _mm256_loadu_si256((const __m256i*)(y + i + 32));
            yl = _mm256_xor_si256(yl, xl);
            yh = _mm256_xor_si256(yh, xh);
            __m256i pl, ph;
            Product16Nibble(yl, yh, lo, hi, pl, ph);
            xl = _mm256_xor_si256(xl, pl);
            xh = _mm256_xor_si256(xh, ph);
            _mm256_storeu_si256((__m256i*)(x + i), xl);
            _mm256_storeu_si256((__m256i*)(x + i + 32), xh);
            _mm256_storeu_si256((__m256i*)(y + i), yl);
            _mm256_storeu_si256((__m256i*)(y + i + 32), yh);
        }
    }
}

__attribute__((target("avx2,gfni")))
inline void Product16Gfni(__m256i dl, __m256i dh, const uint8_t a[4][32],
                          __m256i& pl, __m256i& ph)
{
    pl = _mm256_xor_si256(
        _mm256_gf2p8affine_epi64_epi8(
            dl, _mm256_loadu_si256((const __m256i*)a[0]), 0),
        _mm256_gf2p8affine_epi64_epi8(
            dh, _mm256_loadu_si256((const __m256i*)a[1]), 0));
    ph = _mm256_xor_si256(
        _mm256_gf2p8affine_epi64_epi8(
            dl, _mm256_loadu_si256((const __m256i*)a[2]), 0),
        _mm256_gf2p8affine_epi64_epi8(
            dh, _mm256_loadu_si256((const __m256i*)a[3]), 0));
}

__attribute__((target("avx2,gfni")))
void Gf16SplitGfni(uint8_t* v0, uint8_t* v1, uint8_t* v2, uint8_t* v3,
                   const Map16& m01, const Map16& m23, const Map16& m02,
                   size_t bytes)
{
    struct Pair { uint8_t* x; uint8_t* y; const Map16* m; };
    const Pair layers[4] = {
        { v0, v1, &m01 }, { v2, v3, &m23 }, { v0, v2, &m02 }, { v1, v3, &m02 }
    };
    for (unsigned layer = 0; layer < 4; ++layer)
    {
        const Map16& m = *layers[layer].m;
        uint8_t* x = layers[layer].x;
        uint8_t* y = layers[layer].y;
        for (size_t i = 0; i < bytes; i += 64)
        {
            __m256i xl = _mm256_loadu_si256((const __m256i*)(x + i));
            __m256i xh = _mm256_loadu_si256((const __m256i*)(x + i + 32));
            __m256i yl = _mm256_loadu_si256((const __m256i*)(y + i));
            __m256i yh = _mm256_loadu_si256((const __m256i*)(y + i + 32));
            yl = _mm256_xor_si256(yl, xl);
            yh = _mm256_xor_si256(yh, xh);
            __m256i pl, ph;
            Product16Gfni(yl, yh, m.affine, pl, ph);
            xl = _mm256_xor_si256(xl, pl);
            xh = _mm256_xor_si256(xh, ph);
            _mm256_storeu_si256((__m256i*)(x + i), xl);
            _mm256_storeu_si256((__m256i*)(x + i + 32), xh);
            _mm256_storeu_si256((__m256i*)(y + i), yl);
            _mm256_storeu_si256((__m256i*)(y + i + 32), yh);
        }
    }
}

__attribute__((target("avx2,gfni")))
void Gf16FusedGfni(uint8_t* v0, uint8_t* v1, uint8_t* v2, uint8_t* v3,
                   const Map16& m01, const Map16& m23, const Map16& m02,
                   size_t bytes)
{
    for (size_t i = 0; i < bytes; i += 64)
    {
        __m256i al = _mm256_loadu_si256((const __m256i*)(v0 + i));
        __m256i ah = _mm256_loadu_si256((const __m256i*)(v0 + i + 32));
        __m256i bl = _mm256_loadu_si256((const __m256i*)(v1 + i));
        __m256i bh = _mm256_loadu_si256((const __m256i*)(v1 + i + 32));
        __m256i cl = _mm256_loadu_si256((const __m256i*)(v2 + i));
        __m256i ch = _mm256_loadu_si256((const __m256i*)(v2 + i + 32));
        __m256i dl = _mm256_loadu_si256((const __m256i*)(v3 + i));
        __m256i dh = _mm256_loadu_si256((const __m256i*)(v3 + i + 32));
        __m256i pl, ph;
        bl = _mm256_xor_si256(bl, al);
        bh = _mm256_xor_si256(bh, ah);
        Product16Gfni(bl, bh, m01.affine, pl, ph);
        al = _mm256_xor_si256(al, pl);
        ah = _mm256_xor_si256(ah, ph);
        dl = _mm256_xor_si256(dl, cl);
        dh = _mm256_xor_si256(dh, ch);
        Product16Gfni(dl, dh, m23.affine, pl, ph);
        cl = _mm256_xor_si256(cl, pl);
        ch = _mm256_xor_si256(ch, ph);
        cl = _mm256_xor_si256(cl, al);
        ch = _mm256_xor_si256(ch, ah);
        Product16Gfni(cl, ch, m02.affine, pl, ph);
        al = _mm256_xor_si256(al, pl);
        ah = _mm256_xor_si256(ah, ph);
        dl = _mm256_xor_si256(dl, bl);
        dh = _mm256_xor_si256(dh, bh);
        Product16Gfni(dl, dh, m02.affine, pl, ph);
        bl = _mm256_xor_si256(bl, pl);
        bh = _mm256_xor_si256(bh, ph);
        _mm256_storeu_si256((__m256i*)(v0 + i), al);
        _mm256_storeu_si256((__m256i*)(v0 + i + 32), ah);
        _mm256_storeu_si256((__m256i*)(v1 + i), bl);
        _mm256_storeu_si256((__m256i*)(v1 + i + 32), bh);
        _mm256_storeu_si256((__m256i*)(v2 + i), cl);
        _mm256_storeu_si256((__m256i*)(v2 + i + 32), ch);
        _mm256_storeu_si256((__m256i*)(v3 + i), dl);
        _mm256_storeu_si256((__m256i*)(v3 + i + 32), dh);
    }
}


__attribute__((target("avx2,gfni")))
void Gf8FusedGfniReg(uint8_t* v0, uint8_t* v1, uint8_t* v2, uint8_t* v3,
                     const Map8& m01, const Map8& m23, const Map8& m02,
                     size_t bytes)
{
    const __m256i r01 = _mm256_loadu_si256((const __m256i*)m01.affine);
    const __m256i r23 = _mm256_loadu_si256((const __m256i*)m23.affine);
    const __m256i r02 = _mm256_loadu_si256((const __m256i*)m02.affine);
    for (size_t i = 0; i < bytes; i += 32)
    {
        __m256i a = _mm256_loadu_si256((const __m256i*)(v0 + i));
        __m256i b = _mm256_loadu_si256((const __m256i*)(v1 + i));
        __m256i c = _mm256_loadu_si256((const __m256i*)(v2 + i));
        __m256i d = _mm256_loadu_si256((const __m256i*)(v3 + i));
        b = _mm256_xor_si256(b, a);
        a = _mm256_xor_si256(a, _mm256_gf2p8affine_epi64_epi8(b, r01, 0));
        d = _mm256_xor_si256(d, c);
        c = _mm256_xor_si256(c, _mm256_gf2p8affine_epi64_epi8(d, r23, 0));
        c = _mm256_xor_si256(c, a);
        a = _mm256_xor_si256(a, _mm256_gf2p8affine_epi64_epi8(c, r02, 0));
        d = _mm256_xor_si256(d, b);
        b = _mm256_xor_si256(b, _mm256_gf2p8affine_epi64_epi8(d, r02, 0));
        _mm256_storeu_si256((__m256i*)(v0 + i), a);
        _mm256_storeu_si256((__m256i*)(v1 + i), b);
        _mm256_storeu_si256((__m256i*)(v2 + i), c);
        _mm256_storeu_si256((__m256i*)(v3 + i), d);
    }
}

// GF16 fused hybrid: the twice-used log02 blocks stay in registers, the
// once-used log01/log23 blocks are memory operands.
__attribute__((target("avx2,gfni")))
void Gf16FusedHybrid(uint8_t* v0, uint8_t* v1, uint8_t* v2, uint8_t* v3,
                     const Map16& m01, const Map16& m23, const Map16& m02,
                     size_t bytes)
{
    const __m256i q0 = _mm256_loadu_si256((const __m256i*)m02.affine[0]);
    const __m256i q1 = _mm256_loadu_si256((const __m256i*)m02.affine[1]);
    const __m256i q2 = _mm256_loadu_si256((const __m256i*)m02.affine[2]);
    const __m256i q3 = _mm256_loadu_si256((const __m256i*)m02.affine[3]);
    for (size_t i = 0; i < bytes; i += 64)
    {
        __m256i al = _mm256_loadu_si256((const __m256i*)(v0 + i));
        __m256i ah = _mm256_loadu_si256((const __m256i*)(v0 + i + 32));
        __m256i bl = _mm256_loadu_si256((const __m256i*)(v1 + i));
        __m256i bh = _mm256_loadu_si256((const __m256i*)(v1 + i + 32));
        __m256i cl = _mm256_loadu_si256((const __m256i*)(v2 + i));
        __m256i ch = _mm256_loadu_si256((const __m256i*)(v2 + i + 32));
        __m256i dl = _mm256_loadu_si256((const __m256i*)(v3 + i));
        __m256i dh = _mm256_loadu_si256((const __m256i*)(v3 + i + 32));
        __m256i pl, ph;
        bl = _mm256_xor_si256(bl, al);
        bh = _mm256_xor_si256(bh, ah);
        Product16Gfni(bl, bh, m01.affine, pl, ph);
        al = _mm256_xor_si256(al, pl);
        ah = _mm256_xor_si256(ah, ph);
        dl = _mm256_xor_si256(dl, cl);
        dh = _mm256_xor_si256(dh, ch);
        Product16Gfni(dl, dh, m23.affine, pl, ph);
        cl = _mm256_xor_si256(cl, pl);
        ch = _mm256_xor_si256(ch, ph);
        cl = _mm256_xor_si256(cl, al);
        ch = _mm256_xor_si256(ch, ah);
        pl = _mm256_xor_si256(_mm256_gf2p8affine_epi64_epi8(cl, q0, 0),
                              _mm256_gf2p8affine_epi64_epi8(ch, q1, 0));
        ph = _mm256_xor_si256(_mm256_gf2p8affine_epi64_epi8(cl, q2, 0),
                              _mm256_gf2p8affine_epi64_epi8(ch, q3, 0));
        al = _mm256_xor_si256(al, pl);
        ah = _mm256_xor_si256(ah, ph);
        dl = _mm256_xor_si256(dl, bl);
        dh = _mm256_xor_si256(dh, bh);
        pl = _mm256_xor_si256(_mm256_gf2p8affine_epi64_epi8(dl, q0, 0),
                              _mm256_gf2p8affine_epi64_epi8(dh, q1, 0));
        ph = _mm256_xor_si256(_mm256_gf2p8affine_epi64_epi8(dl, q2, 0),
                              _mm256_gf2p8affine_epi64_epi8(dh, q3, 0));
        bl = _mm256_xor_si256(bl, pl);
        bh = _mm256_xor_si256(bh, ph);
        _mm256_storeu_si256((__m256i*)(v0 + i), al);
        _mm256_storeu_si256((__m256i*)(v0 + i + 32), ah);
        _mm256_storeu_si256((__m256i*)(v1 + i), bl);
        _mm256_storeu_si256((__m256i*)(v1 + i + 32), bh);
        _mm256_storeu_si256((__m256i*)(v2 + i), cl);
        _mm256_storeu_si256((__m256i*)(v2 + i + 32), ch);
        _mm256_storeu_si256((__m256i*)(v3 + i), dl);
        _mm256_storeu_si256((__m256i*)(v3 + i + 32), dh);
    }
}


// GF16 fused radix-four, log02 blocks in registers, log01/log23 blocks
// re-broadcast from the existing 16-byte table rows inside the tile loop.
// Needs no table layout change and no extra memory.
__attribute__((target("avx2,gfni")))
void Gf16FusedBcast(uint8_t* v0, uint8_t* v1, uint8_t* v2, uint8_t* v3,
                    const Map16& m01, const Map16& m23, const Map16& m02,
                    size_t bytes)
{
    const __m256i q0 = _mm256_loadu_si256((const __m256i*)m02.affine[0]);
    const __m256i q1 = _mm256_loadu_si256((const __m256i*)m02.affine[1]);
    const __m256i q2 = _mm256_loadu_si256((const __m256i*)m02.affine[2]);
    const __m256i q3 = _mm256_loadu_si256((const __m256i*)m02.affine[3]);
    for (size_t i = 0; i < bytes; i += 64)
    {
        __m256i al = _mm256_loadu_si256((const __m256i*)(v0 + i));
        __m256i ah = _mm256_loadu_si256((const __m256i*)(v0 + i + 32));
        __m256i bl = _mm256_loadu_si256((const __m256i*)(v1 + i));
        __m256i bh = _mm256_loadu_si256((const __m256i*)(v1 + i + 32));
        __m256i cl = _mm256_loadu_si256((const __m256i*)(v2 + i));
        __m256i ch = _mm256_loadu_si256((const __m256i*)(v2 + i + 32));
        __m256i dl = _mm256_loadu_si256((const __m256i*)(v3 + i));
        __m256i dh = _mm256_loadu_si256((const __m256i*)(v3 + i + 32));
        __m256i pl, ph;
        bl = _mm256_xor_si256(bl, al);
        bh = _mm256_xor_si256(bh, ah);
        pl = _mm256_xor_si256(
            _mm256_gf2p8affine_epi64_epi8(bl, Bcast(m01.affine[0]), 0),
            _mm256_gf2p8affine_epi64_epi8(bh, Bcast(m01.affine[1]), 0));
        ph = _mm256_xor_si256(
            _mm256_gf2p8affine_epi64_epi8(bl, Bcast(m01.affine[2]), 0),
            _mm256_gf2p8affine_epi64_epi8(bh, Bcast(m01.affine[3]), 0));
        al = _mm256_xor_si256(al, pl);
        ah = _mm256_xor_si256(ah, ph);
        dl = _mm256_xor_si256(dl, cl);
        dh = _mm256_xor_si256(dh, ch);
        pl = _mm256_xor_si256(
            _mm256_gf2p8affine_epi64_epi8(dl, Bcast(m23.affine[0]), 0),
            _mm256_gf2p8affine_epi64_epi8(dh, Bcast(m23.affine[1]), 0));
        ph = _mm256_xor_si256(
            _mm256_gf2p8affine_epi64_epi8(dl, Bcast(m23.affine[2]), 0),
            _mm256_gf2p8affine_epi64_epi8(dh, Bcast(m23.affine[3]), 0));
        cl = _mm256_xor_si256(cl, pl);
        ch = _mm256_xor_si256(ch, ph);
        cl = _mm256_xor_si256(cl, al);
        ch = _mm256_xor_si256(ch, ah);
        pl = _mm256_xor_si256(_mm256_gf2p8affine_epi64_epi8(cl, q0, 0),
                              _mm256_gf2p8affine_epi64_epi8(ch, q1, 0));
        ph = _mm256_xor_si256(_mm256_gf2p8affine_epi64_epi8(cl, q2, 0),
                              _mm256_gf2p8affine_epi64_epi8(ch, q3, 0));
        al = _mm256_xor_si256(al, pl);
        ah = _mm256_xor_si256(ah, ph);
        dl = _mm256_xor_si256(dl, bl);
        dh = _mm256_xor_si256(dh, bh);
        pl = _mm256_xor_si256(_mm256_gf2p8affine_epi64_epi8(dl, q0, 0),
                              _mm256_gf2p8affine_epi64_epi8(dh, q1, 0));
        ph = _mm256_xor_si256(_mm256_gf2p8affine_epi64_epi8(dl, q2, 0),
                              _mm256_gf2p8affine_epi64_epi8(dh, q3, 0));
        bl = _mm256_xor_si256(bl, pl);
        bh = _mm256_xor_si256(bh, ph);
        _mm256_storeu_si256((__m256i*)(v0 + i), al);
        _mm256_storeu_si256((__m256i*)(v0 + i + 32), ah);
        _mm256_storeu_si256((__m256i*)(v1 + i), bl);
        _mm256_storeu_si256((__m256i*)(v1 + i + 32), bh);
        _mm256_storeu_si256((__m256i*)(v2 + i), cl);
        _mm256_storeu_si256((__m256i*)(v2 + i + 32), ch);
        _mm256_storeu_si256((__m256i*)(v3 + i), dl);
        _mm256_storeu_si256((__m256i*)(v3 + i + 32), dh);
    }
}

// -------------------------------------------------------------------- harness

struct Group
{
    uint8_t* shard[4];
    size_t bytes;
    ~Group() { for (int i = 0; i < 4; ++i) free(shard[i]); }
};

void Alloc(Group& g, size_t bytes, uint64_t& seed)
{
    g.bytes = bytes;
    for (int i = 0; i < 4; ++i)
    {
        void* p = NULL;
        if (posix_memalign(&p, 64, bytes) != 0) abort();
        uint8_t* d = static_cast<uint8_t*>(p);
        for (size_t j = 0; j < bytes; ++j) d[j] = static_cast<uint8_t>(Rng(seed));
        g.shard[i] = d;
    }
}

template<typename Fn>
double Time(Fn fn, Group& g, unsigned reps)
{
    double best = 1e30;
    for (unsigned r = 0; r < reps; ++r)
    {
        const Clock::time_point t0 = Clock::now();
        fn(g.shard[0], g.shard[1], g.shard[2], g.shard[3], g.bytes);
        const double us =
            std::chrono::duration<double, std::micro>(Clock::now() - t0).count();
        best = std::min(best, us);
    }
    return best;
}

} // namespace

int main(int argc, char** argv)
{
    const unsigned reps = argc > 1 ? unsigned(atoi(argv[1])) : 200;
    uint64_t seed = 0xB5026F5AA96619E9ULL;
    const Map8 a01 = MakeMap8(seed), a23 = MakeMap8(seed), a02 = MakeMap8(seed);
    const Map16 b01 = MakeMap16(seed), b23 = MakeMap16(seed),
                b02 = MakeMap16(seed);

    // equivalence
    {
        const size_t n = 256;
        std::vector<uint8_t> base(n * 4);
        for (size_t i = 0; i < base.size(); ++i)
            base[i] = static_cast<uint8_t>(Rng(seed));
        std::vector<uint8_t> r1(base), r2(base), r3(base);
        uint8_t* p1[4]; uint8_t* p2[4]; uint8_t* p3[4];
        for (int i = 0; i < 4; ++i)
        {
            p1[i] = r1.data() + i * n;
            p2[i] = r2.data() + i * n;
            p3[i] = r3.data() + i * n;
        }
        Gf8SplitNibble(p1[0], p1[1], p1[2], p1[3], a01, a23, a02, n);
        Gf8SplitGfni(p2[0], p2[1], p2[2], p2[3], a01, a23, a02, n);
        Gf8FusedGfni(p3[0], p3[1], p3[2], p3[3], a01, a23, a02, n);
        if (r1 != r2 || r1 != r3) { fprintf(stderr, "GF8 mismatch\n"); return 1; }
        std::vector<uint8_t> r4(base);
        uint8_t* p4[4];
        for (int i = 0; i < 4; ++i) p4[i] = r4.data() + i * n;
        Gf8FusedGfniReg(p4[0], p4[1], p4[2], p4[3], a01, a23, a02, n);
        if (r1 != r4) { fprintf(stderr, "GF8 reg mismatch\n"); return 1; }
        std::vector<uint8_t> r5(base);
        uint8_t* p5[4];
        for (int i = 0; i < 4; ++i) p5[i] = r5.data() + i * n;
        Gf8FusedNibble(p5[0], p5[1], p5[2], p5[3], a01, a23, a02, n);
        if (r1 != r5) { fprintf(stderr, "GF8 fused-nibble mismatch\n"); return 1; }
        r1 = r2 = r3 = base;
        Gf16SplitNibble(p1[0], p1[1], p1[2], p1[3], b01, b23, b02, n);
        Gf16SplitGfni(p2[0], p2[1], p2[2], p2[3], b01, b23, b02, n);
        Gf16FusedGfni(p3[0], p3[1], p3[2], p3[3], b01, b23, b02, n);
        if (r1 != r2 || r1 != r3) { fprintf(stderr, "GF16 mismatch\n"); return 1; }
        r4 = base;
        Gf16FusedHybrid(p4[0], p4[1], p4[2], p4[3], b01, b23, b02, n);
        if (r1 != r4) { fprintf(stderr, "GF16 hybrid mismatch\n"); return 1; }
        r5 = base;
        Gf16FusedBcast(p5[0], p5[1], p5[2], p5[3], b01, b23, b02, n);
        if (r1 != r5) { fprintf(stderr, "GF16 bcast mismatch\n"); return 1; }
    }
    printf("all three schedules agree byte for byte\n\n");

    const size_t sizes[] = { 1024, 4096, 16384, 65536, 262144 };
    printf("GF8 radix-4 inverse group (4 shards), microseconds, best of %u\n",
           reps);
    printf("%10s %12s %12s %12s %12s %10s\n", "bytes", "fused-nibble",
           "gfni-split", "fused-mem", "fused-reg", "gain");
    for (size_t i = 0; i < sizeof(sizes) / sizeof(sizes[0]); ++i)
    {
        Group g; Alloc(g, sizes[i], seed);
        const double s = Time([&](uint8_t* a, uint8_t* b, uint8_t* c,
                                  uint8_t* d, size_t n) {
            Gf8SplitNibble(a, b, c, d, a01, a23, a02, n); }, g, reps);
        const double t = Time([&](uint8_t* a, uint8_t* b, uint8_t* c,
                                  uint8_t* d, size_t n) {
            Gf8SplitGfni(a, b, c, d, a01, a23, a02, n); }, g, reps);
        const double u = Time([&](uint8_t* a, uint8_t* b, uint8_t* c,
                                  uint8_t* d, size_t n) {
            Gf8FusedGfni(a, b, c, d, a01, a23, a02, n); }, g, reps);
        const double w = Time([&](uint8_t* a, uint8_t* b, uint8_t* c,
                                  uint8_t* d, size_t n) {
            Gf8FusedGfniReg(a, b, c, d, a01, a23, a02, n); }, g, reps);
        const double v = Time([&](uint8_t* a, uint8_t* b, uint8_t* c,
                                  uint8_t* d, size_t n) {
            Gf8FusedNibble(a, b, c, d, a01, a23, a02, n); }, g, reps);
        printf("%10zu %12.3f %12.3f %12.3f %12.3f %10.2f (production baseline is fused-nibble)\n",
               sizes[i], v, t, u, w, v / w);
    }

    printf("\nGF16 radix-4 inverse group (4 shards), microseconds\n");
    printf("%10s %12s %12s %12s %12s %10s\n", "bytes", "fused-nibble",
           "gfni-split", "fused-mem", "fused-reg", "gain");
    for (size_t i = 0; i < sizeof(sizes) / sizeof(sizes[0]); ++i)
    {
        Group g; Alloc(g, sizes[i], seed);
        const double s = Time([&](uint8_t* a, uint8_t* b, uint8_t* c,
                                  uint8_t* d, size_t n) {
            Gf16SplitNibble(a, b, c, d, b01, b23, b02, n); }, g, reps);
        const double t = Time([&](uint8_t* a, uint8_t* b, uint8_t* c,
                                  uint8_t* d, size_t n) {
            Gf16SplitGfni(a, b, c, d, b01, b23, b02, n); }, g, reps);
        const double u = Time([&](uint8_t* a, uint8_t* b, uint8_t* c,
                                  uint8_t* d, size_t n) {
            Gf16FusedGfni(a, b, c, d, b01, b23, b02, n); }, g, reps);
        const double w = Time([&](uint8_t* a, uint8_t* b, uint8_t* c,
                                  uint8_t* d, size_t n) {
            Gf16FusedHybrid(a, b, c, d, b01, b23, b02, n); }, g, reps);
        const double v = Time([&](uint8_t* a, uint8_t* b, uint8_t* c,
                                  uint8_t* d, size_t n) {
            Gf16FusedBcast(a, b, c, d, b01, b23, b02, n); }, g, reps);
        printf("%10zu %12.3f %12.3f %12.3f %12.3f %10.2f  bcast=%.3f (%.2fx)\n",
               sizes[i], s, t, u, w, s / w, v, s / v);
    }
    return 0;
}
