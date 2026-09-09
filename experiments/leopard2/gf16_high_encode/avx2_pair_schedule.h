// Experiment only: leopard-79h.38.5.4.18.2. No benchmark clocks.
// Included inside leopard::backend, only by the ordinary AVX2 source variant.
// The caller stores inverse y before consuming its local vectors. No change to
// the table format, transform, scalar tail, or production default is intended.
#ifndef LEO_AVX2_PAIR_SCHEDULE
#define LEO_AVX2_PAIR_SCHEDULE 0
#endif
#if LEO_AVX2_PAIR_SCHEDULE < 0 || LEO_AVX2_PAIR_SCHEDULE > 2
#error Unsupported AVX2 pair schedule
#endif

#if LEO_AVX2_PAIR_SCHEDULE
static LEO_FORCE_INLINE void AVX2FF16PairProductAdd(
    __m256i low_data, __m256i high_data,
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
#if LEO_AVX2_PAIR_SCHEDULE == 2
        // A compiler-only dependency boundary, not a hardware fence. Making
        // subsequent data extraction depend on it limits speculative live
        // nibble/product temporaries. The generated loop must be inspected:
        // reduced spills could still lose to a longer dependency chain.
        __asm__("" : "+x"(low_acc), "+x"(high_acc),
                     "+x"(low_data), "+x"(high_data));
#endif
    }
}
#endif
