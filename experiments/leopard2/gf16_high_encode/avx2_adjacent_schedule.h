// Untimed prototype only: leopard-79h.38.5.4.18.3.
// Included inside leopard::backend by the ordinary AVX2 variant alone.
// Bits: 1=forward pair, 2=accumulating inverse pair. Default is unchanged.
#ifndef LEO_AVX2_ADJACENT_SCHEDULE
#define LEO_AVX2_ADJACENT_SCHEDULE 0
#endif
#if LEO_AVX2_ADJACENT_SCHEDULE < 0 || LEO_AVX2_ADJACENT_SCHEDULE > 3
#error Unsupported adjacent pair schedule
#endif

#if LEO_AVX2_ADJACENT_SCHEDULE
static LEO_FORCE_INLINE void AVX2FF16AdjacentProductAdd(
    __m256i& low_data, __m256i& high_data,
    const __m256i low_tables[4], const __m256i high_tables[4],
    __m256i& low_acc, __m256i& high_acc)
{
    const __m256i mask = _mm256_set1_epi8(15);
    for (unsigned i = 0; i < 4; ++i)
    {
        const __m256i data = i < 2 ? low_data : high_data;
        const __m256i nibble = _mm256_and_si256(
            i % 2 ? _mm256_srli_epi64(data, 4) : data, mask);
        low_acc = _mm256_xor_si256(low_acc,
            _mm256_shuffle_epi8(low_tables[i], nibble));
        high_acc = _mm256_xor_si256(high_acc,
            _mm256_shuffle_epi8(high_tables[i], nibble));
        // Empty GNU asm emits no hardware instruction. Passing data by
        // reference lets the forward caller use these same physical vectors
        // afterward, instead of keeping duplicate pre-boundary y values live.
        __asm__("" : "+x"(low_acc), "+x"(high_acc),
                     "+x"(low_data), "+x"(high_data));
    }
}
#endif
