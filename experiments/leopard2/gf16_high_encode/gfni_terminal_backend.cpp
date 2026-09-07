// Experiment-only .38.5.4.14 member: same GFNI implementation plus one
// accumulating entry. Compile with the original GFNI TU's exact ISA flags.
#include "Leopard2BackendGFNI.cpp"

extern "C" void LeoGFNIFinalAccumulate(
    const void* const* input, void* const* accumulator, unsigned distance,
    uint16_t log01, uint16_t log23, uint16_t log02, uint64_t bytes)
{
    using namespace leopard::backend;
    if (bytes == 0) return;
    for (unsigned i=0; i<distance; ++i)
    {
        // Named pointers avoid SLP-packing a local pointer array using a
        // mnemonic outside the existing GFNI member's reviewed allowlist.
        const auto* value0=static_cast<const uint8_t*>(input[i]);
        const auto* value1=static_cast<const uint8_t*>(input[i+distance]);
        const auto* value2=static_cast<const uint8_t*>(input[i+2U*distance]);
        const auto* value3=static_cast<const uint8_t*>(input[i+3U*distance]);
        auto* sum0=static_cast<uint8_t*>(accumulator[i]);
        auto* sum1=static_cast<uint8_t*>(accumulator[i+distance]);
        auto* sum2=static_cast<uint8_t*>(accumulator[i+2U*distance]);
        auto* sum3=static_cast<uint8_t*>(accumulator[i+3U*distance]);
        uint64_t offset=0;
        while (bytes-offset >= 64)
        {
            __m256i low[4], high[4];
            const auto load=[&](unsigned lane,const uint8_t* value) {
                low[lane]=_mm256_loadu_si256(reinterpret_cast<const __m256i*>(value+offset));
                high[lane]=_mm256_loadu_si256(reinterpret_cast<const __m256i*>(value+offset+32));
            };
            load(0,value0); load(1,value1); load(2,value2); load(3,value3);
            low[1]=_mm256_xor_si256(low[1],low[0]);
            high[1]=_mm256_xor_si256(high[1],high[0]);
            if (log01!=65535)
                AVX2FF16MultiplyAddPair(low[0],high[0],low[1],high[1],FF16Tables[log01]);
            low[3]=_mm256_xor_si256(low[3],low[2]);
            high[3]=_mm256_xor_si256(high[3],high[2]);
            if (log23!=65535)
                AVX2FF16MultiplyAddPair(low[2],high[2],low[3],high[3],FF16Tables[log23]);
            low[2]=_mm256_xor_si256(low[2],low[0]);
            high[2]=_mm256_xor_si256(high[2],high[0]);
            low[3]=_mm256_xor_si256(low[3],low[1]);
            high[3]=_mm256_xor_si256(high[3],high[1]);
            if (log02!=65535)
            {
                AVX2FF16MultiplyAddPair(low[0],high[0],low[2],high[2],FF16Tables[log02]);
                AVX2FF16MultiplyAddPair(low[1],high[1],low[3],high[3],FF16Tables[log02]);
            }
            const auto accumulate=[&](unsigned lane,uint8_t* sum) {
                const __m256i out_low=_mm256_xor_si256(low[lane],
                    _mm256_loadu_si256(reinterpret_cast<const __m256i*>(sum+offset)));
                const __m256i out_high=_mm256_xor_si256(high[lane],
                    _mm256_loadu_si256(reinterpret_cast<const __m256i*>(sum+offset+32)));
                _mm256_storeu_si256(reinterpret_cast<__m256i*>(sum+offset),out_low);
                _mm256_storeu_si256(reinterpret_cast<__m256i*>(sum+offset+32),out_high);
            };
            accumulate(0,sum0); accumulate(1,sum1); accumulate(2,sum2); accumulate(3,sum3);
            offset+=64;
        }
        const uint64_t residual=bytes-offset;
        if (residual)
        {
            // Existing compact GF16 layout/zero-skew handling, no overread.
            alignas(32) uint8_t tail[4][64];
            std::memcpy(tail[0],value0+offset,static_cast<size_t>(residual));
            std::memcpy(tail[1],value1+offset,static_cast<size_t>(residual));
            std::memcpy(tail[2],value2+offset,static_cast<size_t>(residual));
            std::memcpy(tail[3],value3+offset,static_cast<size_t>(residual));
            AVX2FF16Butterfly4Split<true>(tail[0],tail[1],tail[2],tail[3],log01,log23,log02,residual);
            AVX2XorMemory(sum0+offset,tail[0],residual);
            AVX2XorMemory(sum1+offset,tail[1],residual);
            AVX2XorMemory(sum2+offset,tail[2],residual);
            AVX2XorMemory(sum3+offset,tail[3],residual);
        }
    }
}
