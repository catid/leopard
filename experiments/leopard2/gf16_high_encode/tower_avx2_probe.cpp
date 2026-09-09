// Experiment-only pure AVX2 product/conversion; leopard-79h.18.20.2.
#include "tower_avx2_probe.h"
#include <immintrin.h>

static inline __m256i broadcast(const uint8_t* row)
{
    return _mm256_broadcastsi128_si256(_mm_loadu_si128(reinterpret_cast<const __m128i*>(row)));
}

static inline __m256i multiply(__m256i value, __m256i lo, __m256i hi, __m256i mask)
{
    return _mm256_xor_si256(
        _mm256_shuffle_epi8(lo, _mm256_and_si256(value, mask)),
        _mm256_shuffle_epi8(hi, _mm256_and_si256(_mm256_srli_epi16(value, 4), mask)));
}

extern "C" void tower_product_blocks(const uint8_t* source, uint8_t* destination,
                                    size_t blocks, const TowerProductTables* tables)
{
    if (!blocks) return;
    const __m256i mask = _mm256_set1_epi8(15);
    const __m256i t0 = broadcast(tables->row[0]), t1 = broadcast(tables->row[1]);
    const __m256i t2 = broadcast(tables->row[2]), t3 = broadcast(tables->row[3]);
    const __m256i t4 = broadcast(tables->row[4]), t5 = broadcast(tables->row[5]);
    for (size_t block = 0; block < blocks; ++block) {
        const __m256i a = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(source));
        const __m256i b = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(source + 32));
        const __m256i ac = multiply(a, t0, t1, mask);
        const __m256i b_delta_d = multiply(b, t2, t3, mask);
        const __m256i cross = multiply(_mm256_xor_si256(a, b), t4, t5, mask);
        _mm256_storeu_si256(reinterpret_cast<__m256i*>(destination), _mm256_xor_si256(ac, b_delta_d));
        _mm256_storeu_si256(reinterpret_cast<__m256i*>(destination + 32), _mm256_xor_si256(cross, ac));
        source += 64;
        destination += 64;
    }
}

extern "C" void tower_convert_blocks(const uint8_t* source, uint8_t* destination,
                                    size_t blocks, const TowerConvertTables* tables)
{
    if (!blocks) return;
    const __m256i mask = _mm256_set1_epi8(15);
    const __m256i t0 = broadcast(tables->row[0]), t1 = broadcast(tables->row[1]);
    const __m256i t2 = broadcast(tables->row[2]), t3 = broadcast(tables->row[3]);
    for (size_t block = 0; block < blocks; ++block) {
        const __m256i lo = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(source));
        const __m256i hi = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(source + 32));
        const __m256i n0 = _mm256_and_si256(hi, mask);
        const __m256i n1 = _mm256_and_si256(_mm256_srli_epi16(hi, 4), mask);
        const __m256i low_map = _mm256_xor_si256(_mm256_shuffle_epi8(t0, n0), _mm256_shuffle_epi8(t1, n1));
        const __m256i high_map = _mm256_xor_si256(_mm256_shuffle_epi8(t2, n0), _mm256_shuffle_epi8(t3, n1));
        _mm256_storeu_si256(reinterpret_cast<__m256i*>(destination), _mm256_xor_si256(lo, low_map));
        _mm256_storeu_si256(reinterpret_cast<__m256i*>(destination + 32), high_map);
        source += 64;
        destination += 64;
    }
}
