// Experiment-only full pair arithmetic, leopard-79h.18.20.3.
#include "tower_butterfly_probe.h"
#include <immintrin.h>

namespace {
inline __m256i broadcast(const uint8_t* row)
{
    return _mm256_broadcastsi128_si256(_mm_loadu_si128(reinterpret_cast<const __m128i*>(row)));
}
inline __m256i load(const uint8_t* p) { return _mm256_loadu_si256(reinterpret_cast<const __m256i*>(p)); }
inline void store(uint8_t* p, __m256i x) { _mm256_storeu_si256(reinterpret_cast<__m256i*>(p), x); }
inline uint8_t small(const uint8_t* lo, const uint8_t* hi, uint8_t value)
{
    return lo[value & 15] ^ hi[value >> 4];
}
inline uint16_t scalar_product(uint16_t value, const TowerProductTables& t)
{
    const uint8_t a = uint8_t(value), b = uint8_t(value >> 8);
    const uint8_t ac = small(t.row[0], t.row[1], a);
    const uint8_t bd = small(t.row[2], t.row[3], b);
    const uint8_t cross = small(t.row[4], t.row[5], a ^ b);
    return uint16_t(ac ^ bd) | (uint16_t(cross ^ ac) << 8);
}
inline __m256i small_vector(__m256i x, __m256i low, __m256i high, __m256i mask)
{
    return _mm256_xor_si256(
        _mm256_shuffle_epi8(low, _mm256_and_si256(x, mask)),
        _mm256_shuffle_epi8(high, _mm256_and_si256(_mm256_srli_epi16(x, 4), mask)));
}

template<bool Inverse, bool Add, bool Zero>
void pair(const uint8_t* x, const uint8_t* y, uint8_t* u, uint8_t* v,
          size_t bytes, const TowerProductTables* tables)
{
    // Dispatch checks bytes before entry, so even null tables are valid for
    // zero-byte and zero-skew operations. Zero is a compile-time branch.
    const __m256i mask = _mm256_set1_epi8(15);
    __m256i t0{}, t1{}, t2{}, t3{}, t4{}, t5{};
    if (!Zero) {
        t0 = broadcast(tables->row[0]); t1 = broadcast(tables->row[1]);
        t2 = broadcast(tables->row[2]); t3 = broadcast(tables->row[3]);
        t4 = broadcast(tables->row[4]); t5 = broadcast(tables->row[5]);
    }
    size_t offset = 0;
    for (; bytes - offset >= 64; offset += 64) {
        __m256i xl = load(x + offset), xh = load(x + offset + 32);
        __m256i yl = load(y + offset), yh = load(y + offset + 32);
        if (Inverse) {
            yl = _mm256_xor_si256(yl, xl); yh = _mm256_xor_si256(yh, xh);
        }
        if (!Zero) {
            const __m256i ac = small_vector(yl, t0, t1, mask);
            const __m256i bd = small_vector(yh, t2, t3, mask);
            const __m256i cross = small_vector(_mm256_xor_si256(yl, yh), t4, t5, mask);
            xl = _mm256_xor_si256(xl, _mm256_xor_si256(ac, bd));
            xh = _mm256_xor_si256(xh, _mm256_xor_si256(cross, ac));
        }
        if (!Inverse) {
            yl = _mm256_xor_si256(yl, xl); yh = _mm256_xor_si256(yh, xh);
        }
        if (Add) {
            xl = _mm256_xor_si256(xl, load(u + offset));
            xh = _mm256_xor_si256(xh, load(u + offset + 32));
            yl = _mm256_xor_si256(yl, load(v + offset));
            yh = _mm256_xor_si256(yh, load(v + offset + 32));
        }
        store(u + offset, xl); store(u + offset + 32, xh);
        store(v + offset, yl); store(v + offset + 32, yh);
    }
    const size_t symbols = (bytes - offset) / 2;
    for (size_t i = 0; i < symbols; ++i) {
        uint16_t a = x[offset+i] | (uint16_t(x[offset+symbols+i]) << 8);
        uint16_t b = y[offset+i] | (uint16_t(y[offset+symbols+i]) << 8);
        if (Inverse) b ^= a;
        if (!Zero) a ^= scalar_product(b, *tables);
        if (!Inverse) b ^= a;
        if (Add) {
            a ^= u[offset+i] | (uint16_t(u[offset+symbols+i]) << 8);
            b ^= v[offset+i] | (uint16_t(v[offset+symbols+i]) << 8);
        }
        u[offset+i] = uint8_t(a); u[offset+symbols+i] = uint8_t(a >> 8);
        v[offset+i] = uint8_t(b); v[offset+symbols+i] = uint8_t(b >> 8);
    }
}

template<bool Inverse, bool Add>
void dispatch(const uint8_t* x, const uint8_t* y, uint8_t* u, uint8_t* v,
              size_t bytes, uint16_t log, const TowerProductTables* tables)
{
    if (!bytes) return;
    if (log == 65535) pair<Inverse, Add, true>(x, y, u, v, bytes, nullptr);
    else pair<Inverse, Add, false>(x, y, u, v, bytes, tables);
}
} // namespace

extern "C" void tower_convert_involution(const uint8_t* input, uint8_t* output,
                                        size_t bytes, const TowerLowMapTables* map)
{
    if (!bytes) return;
    const __m256i t0 = broadcast(map->row[0]), t1 = broadcast(map->row[1]);
    const __m256i mask = _mm256_set1_epi8(15);
    size_t offset = 0;
    for (; bytes - offset >= 64; offset += 64) {
        const __m256i lo = load(input + offset), hi = load(input + offset + 32);
        const __m256i mapped = small_vector(hi, t0, t1, mask);
        store(output + offset, _mm256_xor_si256(lo, mapped));
        store(output + offset + 32, hi);
    }
    const size_t symbols = (bytes - offset) / 2;
    for (size_t i = 0; i < symbols; ++i) {
        const uint8_t hi = input[offset+symbols+i];
        const uint8_t lo = input[offset+i] ^ small(map->row[0], map->row[1], hi);
        output[offset+i] = lo; output[offset+symbols+i] = hi;
    }
}
extern "C" void tower_ifft_pair(uint8_t* x, uint8_t* y, size_t bytes,
                               uint16_t log, const TowerProductTables* t)
{
    dispatch<true, false>(x, y, x, y, bytes, log, t);
}
extern "C" void tower_fft_pair(uint8_t* x, uint8_t* y, size_t bytes,
                              uint16_t log, const TowerProductTables* t)
{
    dispatch<false, false>(x, y, x, y, bytes, log, t);
}
extern "C" void tower_fft_out(const uint8_t* x, const uint8_t* y,
                             uint8_t* u, uint8_t* v, size_t bytes,
                             uint16_t log, const TowerProductTables* t)
{
    dispatch<false, false>(x, y, u, v, bytes, log, t);
}
extern "C" void tower_ifft_accumulate(const uint8_t* x, const uint8_t* y,
                                     uint8_t* u, uint8_t* v, size_t bytes,
                                     uint16_t log, const TowerProductTables* t)
{
    dispatch<true, true>(x, y, u, v, bytes, log, t);
}
